"""MarkerTool — click on a body to attach a marker at the click point."""

from __future__ import annotations

import math

from PySide6.QtCore import Qt

from ..models import MarkerSpec
from ..widgets import BodyItem
from .tool_base import Tool


class MarkerTool(Tool):
    """Click on a body → adds a marker at the clicked local point.

    If the click misses every body, the marker is attached to the
    implicit ``ground`` body at the (global) click location.
    """

    name = "marker"
    cursor = Qt.CrossCursor

    def activate(self) -> None:
        super().activate()
        self.window.statusBar().showMessage(
            "Marker: click on a body to attach a marker (click on empty "
            "space → attach to ground; Esc to cancel)…", 0)

    def deactivate(self) -> None:
        super().deactivate()
        self.window.statusBar().clearMessage()

    def mouse_press(self, event) -> bool:
        if event.button() != Qt.LeftButton:
            return False
        scene_pt = event.scenePos()
        snapped  = self.scene.snap(scene_pt)

        # Find topmost body item under the cursor.
        body_item = None
        for it in self.scene.items(scene_pt):
            if isinstance(it, BodyItem):
                body_item = it
                break

        if body_item is not None:
            owning_body = body_item.spec
            local = body_item.mapFromScene(snapped[0], snapped[1])
            local_xy = (float(local.x()), float(local.y()))
            # The marker inherits the parent body's orientation so that
            # its axes line up with the body frame at creation time.
            theta = float(owning_body.orientation)
        else:
            owning_body = next(
                (b for b in self.spec.bodies if b.id == "ground"), None)
            local_xy = (snapped[0], snapped[1])
            theta = 0.0

        body_id   = owning_body.id   if owning_body is not None else "ground"
        body_name = owning_body.name if owning_body is not None else "ground"
        n_on_body = sum(1 for m in self.spec.markers if m.body_id == body_id) + 1
        spec = MarkerSpec(
            name=f"{body_name}.mk{n_on_body}",
            body_id=body_id,
            local_position=local_xy,
            theta=theta,
        )
        self.spec.markers.append(spec)
        self.window.add_marker_item(spec)
        self._commit()
        # One-shot: revert to Select after placing the marker.
        self.window.set_active_tool("select")
        return True
