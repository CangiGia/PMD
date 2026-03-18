"""
Planar Multi-Body Dynamics Topological Assembler

This module implements the _assemble() function which computes body
positions and orientations from joint topology and q0 values.

Author: Giacomo Cangi
"""

import numpy as np
from collections import deque
from .model import Ground, Marker
from .mechanics import A_matrix, s_rot


def _assemble(bodies, joints):
    """Compute r and p for all bodies from joint topology and q0 values."""

    # Build spanning tree via BFS from Ground
    visited = {id(Ground): Ground}
    queue = deque([Ground])
    spanning_tree = []   # list of (parent_body, child_body, joint)
    loop_joints = []
    spanning_set = set() # set of joint ids already in spanning_tree

    while queue:
        node = queue.popleft()
        for joint in joints:
            if joint.iMarker is None or joint.jMarker is None:
                continue   # disc, rel-rot, rel-tran, rigid by bodies: skip
            bi = joint.iMarker.body
            bj = joint.jMarker.body
            jid = id(joint)
            if id(bi) == id(node) and id(bj) not in visited:
                spanning_tree.append((bi, bj, joint))
                spanning_set.add(jid)
                visited[id(bj)] = bj
                queue.append(bj)
            elif id(bj) == id(node) and id(bi) not in visited:
                spanning_tree.append((bj, bi, joint))
                spanning_set.add(jid)
                visited[id(bi)] = bi
                queue.append(bi)
            elif jid not in spanning_set and (
                    (id(bi) in visited and id(bj) in visited)):
                if jid not in [id(lj) for lj in loop_joints]:
                    loop_joints.append(joint)

    # Forward kinematics along spanning tree
    for (parent, child, joint) in spanning_tree:
        _fk_step(parent, child, joint)

    # Resolve deferred markers
    _resolve_deferred_markers(bodies)

    # Newton-Raphson loop closing
    if loop_joints:
        _close_loops(loop_joints, bodies)


def _fk_step(parent, child, joint):
    """Place child body given parent state and joint.q0."""
    q0 = joint.q0
    im = joint.iMarker
    jm = joint.jMarker
    # identify which marker belongs to parent and which to child
    if id(im.body) == id(parent):
        parent_marker, child_marker = im, jm
    else:
        parent_marker, child_marker = jm, im

    parent._A = A_matrix(parent.p)

    if joint.type == 'rev':
        # global joint position from parent side
        parent_sP = parent._A @ parent_marker.position.reshape(2, 1)
        r_joint = parent.r + parent_sP
        # child orientation
        child.p = parent.p + q0
        child._A = A_matrix(child.p)
        # child position
        child_sP = child._A @ child_marker.position.reshape(2, 1)
        child.r = r_joint - child_sP

    elif joint.type == 'tran':
        # child keeps same orientation as parent
        child.p = parent.p
        child._A = A_matrix(child.p)
        # slide direction in global frame (use parent marker orientation)
        if parent_marker.has_orientation:
            slide_dir = parent._A @ parent_marker._ulocal
        else:
            slide_dir = child._A @ child_marker._ulocal
        parent_sP = parent._A @ parent_marker.position.reshape(2, 1)
        r_joint = parent.r + parent_sP
        child_sP = child._A @ child_marker.position.reshape(2, 1)
        child.r = r_joint - child_sP + q0 * slide_dir

    elif joint.type == 'rev-rev':
        # q0 = global angle of segment from parent pivot to child pivot
        parent_sP = parent._A @ parent_marker.position.reshape(2, 1)
        r_joint = parent.r + parent_sP
        child._A = A_matrix(child.p)
        child_sP = child._A @ child_marker.position.reshape(2, 1)
        L = joint.L
        r_child_pivot = r_joint + L * np.array([[np.cos(q0)], [np.sin(q0)]])
        child.r = r_child_pivot - child_sP

    # disc, rel-rot, rel-tran, rigid: no positional FK needed


def _resolve_deferred_markers(bodies):
    """Compute local position of markers created with add_marker_at()."""
    for body in bodies:
        for m in body._markers:
            if not hasattr(m, '_deferred_ref'):
                continue
            ref = m._deferred_ref
            offset = m._deferred_offset
            ref_body = ref.body
            if ref_body is Ground:
                ref_rP = ref.position.reshape(2, 1)
                ref_A = np.eye(2)
            else:
                ref_body._A = A_matrix(ref_body.p)
                ref_sP = ref_body._A @ ref.position.reshape(2, 1)
                ref_rP = ref_body.r + ref_sP
                ref_A = ref_body._A
            # global position = ref_rP + R(ref_body) @ offset
            offset_global = ref_A @ np.asarray(offset, dtype=float).reshape(2, 1)
            new_rP_global = ref_rP + offset_global
            # convert to local frame of m.body
            child_body = m.body
            child_body._A = A_matrix(child_body.p)
            local_col = child_body._A.T @ (new_rP_global - child_body.r)
            m.position = local_col.flatten()
            del m._deferred_ref
            del m._deferred_offset


def _eval_loop_phi(loop_joints):
    rows = []
    for joint in loop_joints:
        im = joint.iMarker
        jm = joint.jMarker
        bi = im.body
        bj = jm.body
        bi_A = np.eye(2) if bi is Ground else A_matrix(bi.p)
        bj_A = np.eye(2) if bj is Ground else A_matrix(bj.p)
        ri = (bi.r if bi is not Ground else np.zeros((2,1))) + bi_A @ im.position.reshape(2,1)
        rj = (bj.r if bj is not Ground else np.zeros((2,1))) + bj_A @ jm.position.reshape(2,1)
        if joint.type == 'rev':
            rows.extend((ri - rj).flatten().tolist())
    return np.array(rows)


def _eval_loop_jacobian(loop_joints, bodies):
    nq = 3 * len(bodies)
    body_idx = {id(b): i for i, b in enumerate(bodies)}
    rows = []
    for joint in loop_joints:
        if joint.type != 'rev':
            continue
        im = joint.iMarker
        jm = joint.jMarker
        bi = im.body
        bj = jm.body
        row_x = np.zeros(nq)
        row_y = np.zeros(nq)
        if bi is not Ground and id(bi) in body_idx:
            ci = body_idx[id(bi)]
            bi_A = A_matrix(bi.p)
            si = bi_A @ im.position.reshape(2,1)
            sir = s_rot(si)
            row_x[3*ci]   += 1;  row_x[3*ci+2] += sir[0,0]
            row_y[3*ci+1] += 1;  row_y[3*ci+2] += sir[1,0]
        if bj is not Ground and id(bj) in body_idx:
            cj = body_idx[id(bj)]
            bj_A = A_matrix(bj.p)
            sj = bj_A @ jm.position.reshape(2,1)
            sjr = s_rot(sj)
            row_x[3*cj]   -= 1;  row_x[3*cj+2] -= sjr[0,0]
            row_y[3*cj+1] -= 1;  row_y[3*cj+2] -= sjr[1,0]
        rows.append(row_x)
        rows.append(row_y)
    return np.array(rows)


def _close_loops(loop_joints, bodies):
    """Newton-Raphson to satisfy loop-closing constraints."""
    for iteration in range(50):
        phi = _eval_loop_phi(loop_joints)
        if np.linalg.norm(phi) < 1e-12:
            return
        J = _eval_loop_jacobian(loop_joints, bodies)
        dq, _, _, _ = np.linalg.lstsq(J, -phi, rcond=None)
        for Bi, body in enumerate(bodies):
            body.r = body.r + dq[3*Bi:3*Bi+2].reshape(2, 1)
            body.p = body.p + float(dq[3*Bi+2])
            body._A = A_matrix(body.p)
        _resolve_deferred_markers(bodies)
    raise RuntimeError(
        "Assembly Newton-Raphson did not converge. "
        "Check model geometry and q0 values.")
