"""Interactive editing tools for the pre-processor canvas."""

from .tool_base import Tool
from .select_tool import SelectTool
from .body_rect_tool import BodyRectTool
from .body_link_tool import BodyLinkTool
from .marker_tool import MarkerTool
from .joint_tool import JointTool, RevJointTool, TranJointTool
from .ptp_force_tool import PtpForceTool


# Registry: tool key (matches RibbonBar) → factory.
TOOL_REGISTRY = {
    "select":      SelectTool,
    "body_rect":   BodyRectTool,
    "body_link":   BodyLinkTool,
    "marker":      MarkerTool,
    "joint_rev":   RevJointTool,
    "joint_tran":  TranJointTool,
    "force_ptp":   PtpForceTool,
}


def make_tool(key: str, window) -> Tool:
    """Instantiate the tool registered under ``key`` (fallback: SelectTool)."""
    cls = TOOL_REGISTRY.get(key, SelectTool)
    return cls(window)


__all__ = [
    "Tool", "SelectTool",
    "BodyRectTool", "BodyLinkTool",
    "MarkerTool",
    "JointTool", "RevJointTool", "TranJointTool",
    "PtpForceTool",
    "TOOL_REGISTRY", "make_tool",
]
