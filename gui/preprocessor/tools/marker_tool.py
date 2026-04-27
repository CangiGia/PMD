"""MarkerTool — click on a body to attach a marker at the click point."""

from __future__ import annotations

import math

from PySide6.QtCore import Qt

from ..models import MarkerSpec
from ..widgets import BodyItem
from .tool_base import Tool


class MarkerTool(Tool):
    """Click on a body → adds a marker at the clicked local point.

    If the click is on empty space, the marker is attached to ``ground``
    at the (global) click location.
    """

    name = "marker"
    cursor = Qt.CrossCursor

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
            local = body_item.mapFromScene(snapped[0], snapped[1])
            spec = MarkerSpec(
                name=f"mrk{len(self.spec.markers) + 1}",
                body_id=body_item.spec.id,
                local_position=(float(local.x()), float(local.y())),
                theta=0.0,
            )
        else:
            spec = MarkerSpec(
                name=f"mrk{len(self.spec.markers) + 1}",
                body_id="ground",
                local_position=(snapped[0], snapped[1]),
                theta=0.0,
            )
        self.spec.markers.append(spec)
        self.window.add_marker_item(spec)
        self._commit()
        return True
