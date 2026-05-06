"""Session — container for a single solved pmd simulation."""

import inspect
import pathlib

_PMD_DIR = pathlib.Path(__file__).parent.parent.resolve()


class Session:
    """Wraps a solved PlanarMultibodyModel together with its time results.

    Parameters
    ----------
    model : PlanarMultibodyModel
        The solved model (_distribute_results already called).
    T : ndarray
        Time vector, shape (nSteps,).
    uT : ndarray
        State matrix, shape (nSteps, 2*nB3).
    name : str, optional
        Human-readable label shown in the UI.
        Defaults to "Simulation N" where N is an auto-incremented counter.
    """

    _count = 0

    def __init__(self, model, T, uT, name=None, preprocessor_spec=None):
        Session._count += 1
        self.model = model
        self.T = T
        self.uT = uT
        # Optional: when the simulation was launched from the
        # preprocessor we keep the originating ModelSpec so the
        # animation pane can reuse the preprocessor's body / marker /
        # joint visuals instead of falling back to coloured circles.
        self.preprocessor_spec = preprocessor_spec
        if name is None:
            try:
                for frame_info in inspect.stack()[1:]:
                    caller = pathlib.Path(frame_info.filename).resolve()
                    if not caller.is_relative_to(_PMD_DIR):
                        name = caller.stem
                        break
            except Exception:
                pass
            if name is None:
                name = f"Simulation {Session._count}"
        self.name = name

    def __repr__(self):
        n_bodies = len(self.model.Bodies)
        n_joints = len(self.model.Joints)
        n_steps = len(self.T)
        return (
            f"Session('{self.name}', "
            f"bodies={n_bodies}, joints={n_joints}, steps={n_steps})"
        )
