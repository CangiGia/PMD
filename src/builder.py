"""Planar Multi-Body Dynamics Topological Assembler.

This module implements the ``_assemble()`` function which computes body
positions and orientations from joint topology and q0 values.

Author: Giacomo Cangi
"""

import numpy as np
from collections import deque

from .model import Ground, Marker
from .mechanics import rotation_matrix, rotate_90


def _assemble(bodies, joints):
    """Compute position and orientation for all bodies from joint topology.

    Uses BFS spanning tree from Ground, forward kinematics along tree
    edges, deferred marker resolution, and Newton-Raphson loop closing.

    Parameters
    ----------
    bodies : list of Body
        List of Body objects.
    joints : list of Joint
        List of Joint objects.
    """
    # Build spanning tree via BFS from Ground
    visited = {id(Ground): Ground}
    queue = deque([Ground])
    spanning_tree = []
    loop_joints = []
    spanning_set = set()

    while queue:
        node = queue.popleft()
        for joint in joints:
            if joint.iMarker is None or joint.jMarker is None:
                continue
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

    # Forward kinematics along spanning tree — polymorphic dispatch
    for (parent, child, joint) in spanning_tree:
        joint.fk_step(parent, child)

    # Resolve deferred markers
    _resolve_deferred_markers(bodies)

    # Newton-Raphson loop closing
    if loop_joints:
        _close_loops(loop_joints, bodies)


def _resolve_deferred_markers(bodies):
    """Compute local position of markers created with ``add_marker_at()``.

    Parameters
    ----------
    bodies : list of Body
        List of Body objects.
    """
    for body in bodies:
        for m in body._markers:
            if not hasattr(m, '_deferred_ref'):
                continue
            ref = m._deferred_ref
            offset = m._deferred_offset
            ref_body = ref.body
            if ref_body is Ground:
                ref_rP = ref.local_position.reshape(2, 1)
                ref_A = np.eye(2)
            else:
                ref_body._rotation_matrix = rotation_matrix(ref_body.orientation)
                ref_sP = ref_body._rotation_matrix @ ref.local_position.reshape(2, 1)
                ref_rP = ref_body.position + ref_sP
                ref_A = ref_body._rotation_matrix
            offset_global = ref_A @ np.asarray(offset, dtype=float).reshape(2, 1)
            new_rP_global = ref_rP + offset_global
            child_body = m.body
            child_body._rotation_matrix = rotation_matrix(child_body.orientation)
            local_col = child_body._rotation_matrix.T @ (new_rP_global - child_body.position)
            m.local_position = local_col.flatten()
            del m._deferred_ref
            del m._deferred_offset


def _eval_loop_phi(loop_joints):
    """Evaluate loop-closing constraint residuals.

    Parameters
    ----------
    loop_joints : list of Joint
        List of loop-closing Joint objects.

    Returns
    -------
    ndarray
        1-D array of constraint residuals.
    """
    rows = []
    for joint in loop_joints:
        im = joint.iMarker
        jm = joint.jMarker
        bi = im.body
        bj = jm.body
        bi_A = np.eye(2) if bi is Ground else rotation_matrix(bi.orientation)
        bj_A = np.eye(2) if bj is Ground else rotation_matrix(bj.orientation)
        ri = ((bi.position if bi is not Ground else np.zeros((2, 1)))
              + bi_A @ im.local_position.reshape(2, 1))
        rj = ((bj.position if bj is not Ground else np.zeros((2, 1)))
              + bj_A @ jm.local_position.reshape(2, 1))
        from .constraints import RevJoint
        if isinstance(joint, RevJoint):
            rows.extend((ri - rj).flatten().tolist())
    return np.array(rows)


def _eval_loop_jacobian(loop_joints, bodies):
    """Evaluate loop-closing constraint Jacobian.

    Parameters
    ----------
    loop_joints : list of Joint
        List of loop-closing Joint objects.
    bodies : list of Body
        List of Body objects.

    Returns
    -------
    ndarray
        2-D array of shape (nRows, 3*nBodies).
    """
    nq = 3 * len(bodies)
    body_idx = {id(b): i for i, b in enumerate(bodies)}
    rows = []
    from .constraints import RevJoint
    for joint in loop_joints:
        if not isinstance(joint, RevJoint):
            continue
        im = joint.iMarker
        jm = joint.jMarker
        bi = im.body
        bj = jm.body
        row_x = np.zeros(nq)
        row_y = np.zeros(nq)
        if bi is not Ground and id(bi) in body_idx:
            ci = body_idx[id(bi)]
            bi_A = rotation_matrix(bi.orientation)
            si = bi_A @ im.local_position.reshape(2, 1)
            sir = rotate_90(si)
            row_x[3*ci] += 1
            row_x[3*ci+2] += sir[0, 0]
            row_y[3*ci+1] += 1
            row_y[3*ci+2] += sir[1, 0]
        if bj is not Ground and id(bj) in body_idx:
            cj = body_idx[id(bj)]
            bj_A = rotation_matrix(bj.orientation)
            sj = bj_A @ jm.local_position.reshape(2, 1)
            sjr = rotate_90(sj)
            row_x[3*cj] -= 1
            row_x[3*cj+2] -= sjr[0, 0]
            row_y[3*cj+1] -= 1
            row_y[3*cj+2] -= sjr[1, 0]
        rows.append(row_x)
        rows.append(row_y)
    return np.array(rows)


def _close_loops(loop_joints, bodies):
    """Newton-Raphson to satisfy loop-closing constraints.

    Parameters
    ----------
    loop_joints : list of Joint
        List of loop-closing Joint objects.
    bodies : list of Body
        List of Body objects.

    Raises
    ------
    RuntimeError
        If convergence is not achieved in 50 iterations.
    """
    for iteration in range(50):
        phi = _eval_loop_phi(loop_joints)
        if np.linalg.norm(phi) < 1e-12:
            return
        J = _eval_loop_jacobian(loop_joints, bodies)
        dq, _, _, _ = np.linalg.lstsq(J, -phi, rcond=None)
        for Bi, body in enumerate(bodies):
            body.position = body.position + dq[3*Bi:3*Bi+2].reshape(2, 1)
            body.orientation = body.orientation + float(dq[3*Bi+2])
            body._rotation_matrix = rotation_matrix(body.orientation)
        _resolve_deferred_markers(bodies)
    raise RuntimeError(
        "Assembly Newton-Raphson did not converge. "
        "Check model geometry and q0 values.")
