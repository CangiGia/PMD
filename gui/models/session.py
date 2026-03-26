"""Session — container for a single solved PMD simulation."""


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

    def __init__(self, model, T, uT, name=None):
        Session._count += 1
        self.model = model
        self.T = T
        self.uT = uT
        self.name = name or f"Simulation {Session._count}"

    def __repr__(self):
        n_bodies = len(self.model.Bodies)
        n_joints = len(self.model.Joints)
        n_steps = len(self.T)
        return (
            f"Session('{self.name}', "
            f"bodies={n_bodies}, joints={n_joints}, steps={n_steps})"
        )
