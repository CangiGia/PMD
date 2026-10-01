# GUI

PMD includes a GUI tool built with PyQt6.

## PostProcessor

Visualise simulation results after calling `model.solve()`:

```python
from pmd.gui import PostProcessor
pp = PostProcessor(model)
pp.show()
```

The post-processor provides:

- Animated trajectory playback.
- Time-history plots for position, velocity, and acceleration of each body.
- Joint reaction force plots.
- Export to CSV or image.

## Theme helpers

```python
from pmd.gui import apply_light_theme, apply_dark_theme

apply_dark_theme(app)   # pass a QApplication instance
```
