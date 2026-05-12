# GUI

PMD includes two GUI tools built with PyQt6.

## PreProcessor

Interactive model builder.  Launch with:

```python
from pmd.gui import PreProcessor
PreProcessor().show()
```

The pre-processor lets you:

- Add and position bodies graphically.
- Attach markers and joints via point-and-click.
- Configure force elements and driver functions.
- Export the assembled model as a JSON file.

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
