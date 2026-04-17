import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerPatch
import matplotlib.transforms

# custom legend elements

# styling for smaller figures (linewidths are reduced)

ADJUST_X = 0.7
ADJUST_Y = 0.7
LINEWIDTH = 0.5
SLASHWIDTH = 0.75

class ForwardSlashRectangle(Rectangle):
    def __init__(self, xy=(0,0), width=0, height=0, *, angle=0.0, rotation_point='xy', **kwargs):
        super().__init__(xy, width, height, angle=angle, rotation_point=rotation_point, **kwargs)
        
class BackSlashRectangle(Rectangle):
    def __init__(self, xy=(0,0), width=0, height=0, *, angle=0.0, rotation_point='xy', **kwargs):
        super().__init__(xy, width, height, angle=angle, rotation_point=rotation_point, **kwargs)
        
class XRectangle(Rectangle):
    def __init__(self, xy=(0,0), width=0, height=0, *, angle=0.0, rotation_point='xy', **kwargs):
        super().__init__(xy, width, height, angle=angle, rotation_point=rotation_point, **kwargs)

class DotRectangle(Rectangle):
    def __init__(self, xy=(0,0), width=0, height=0, *, angle=0.0, rotation_point='xy', **kwargs):
        super().__init__(xy, width, height, angle=angle, rotation_point=rotation_point, **kwargs)
    
class BackSlashHandler(HandlerPatch):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        lower_left = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        p = Rectangle(xy=lower_left, width=width, height=height, facecolor='white', edgecolor='black', linewidth=LINEWIDTH)
        handlebox.add_artist(p)
        l = Line2D(xdata=[0 + ADJUST_X, width - handlebox.xdescent], ydata=[0 + ADJUST_Y, height - handlebox.ydescent + ADJUST_Y], linewidth=SLASHWIDTH, color='black')
        handlebox.add_artist(l)
        return p
    
class ForwardSlashHandler(HandlerPatch):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        lower_left = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        p = Rectangle(xy=lower_left, width=width, height=height, facecolor='white', edgecolor='black', linewidth=LINEWIDTH)
        handlebox.add_artist(p)
        l = Line2D(xdata=[0 + ADJUST_X, width - handlebox.xdescent], ydata=[height - handlebox.ydescent + ADJUST_Y, 0 + ADJUST_Y], linewidth=SLASHWIDTH, color='black')
        handlebox.add_artist(l)
        return p

class XHandler(HandlerPatch):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        lower_left = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        p = Rectangle(xy=lower_left, width=width, height=height, facecolor='white', edgecolor='black', linewidth=LINEWIDTH)
        handlebox.add_artist(p)
        l1 = Line2D(xdata=[0 + ADJUST_X, width - handlebox.xdescent], ydata=[height - handlebox.ydescent + ADJUST_Y, 0 + ADJUST_Y], linewidth=SLASHWIDTH, color='black')
        handlebox.add_artist(l1)
        l2 = Line2D(xdata=[0 + ADJUST_X, width - handlebox.xdescent], ydata=[0 + ADJUST_Y, height - handlebox.ydescent + ADJUST_Y], linewidth=SLASHWIDTH, color='black')
        handlebox.add_artist(l2)
        return p
    
class DotRectangleHandler(HandlerPatch):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        lower_left = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        center = 0.5 * (width - handlebox.xdescent + ADJUST_X), 0.5 * (height - handlebox.ydescent + ADJUST_Y)
        p = Rectangle(xy=lower_left, width=width, height=height, facecolor='white', edgecolor='black', linewidth=LINEWIDTH)
        handlebox.add_artist(p)
        c = Circle(xy=center, radius=0.08 * width, facecolor='black')
        handlebox.add_artist(c)
        return p
    
# styling for large figures (linewidths are larger)

LARGE_ADJUST_X = 1
LARGE_ADJUST_Y = 1
LARGE_LINEWIDTH = 1
LARGE_SLASHWIDTH = 1
    
class LargeBackSlashHandler(HandlerPatch):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        lower_left = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        p = Rectangle(xy=lower_left, width=width, height=height, facecolor='white', edgecolor='black', linewidth=LARGE_LINEWIDTH)
        handlebox.add_artist(p)
        l = Line2D(xdata=[0 + LARGE_ADJUST_X, width - handlebox.xdescent], ydata=[0 + LARGE_ADJUST_Y, height - handlebox.ydescent + LARGE_ADJUST_Y], linewidth=LARGE_SLASHWIDTH, color='black')
        handlebox.add_artist(l)
        return p
    
class LargeForwardSlashHandler(HandlerPatch):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        lower_left = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        p = Rectangle(xy=lower_left, width=width, height=height, facecolor='white', edgecolor='black', linewidth=LARGE_LINEWIDTH)
        handlebox.add_artist(p)
        l = Line2D(xdata=[0 + LARGE_ADJUST_X, width - handlebox.xdescent], ydata=[height - handlebox.ydescent + LARGE_ADJUST_Y, 0 + LARGE_ADJUST_Y], linewidth=LARGE_SLASHWIDTH, color='black')
        handlebox.add_artist(l)
        return p

class LargeXHandler(HandlerPatch):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        lower_left = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        p = Rectangle(xy=lower_left, width=width, height=height, facecolor='white', edgecolor='black', linewidth=LARGE_LINEWIDTH)
        handlebox.add_artist(p)
        l1 = Line2D(xdata=[0 + LARGE_ADJUST_X, width - handlebox.xdescent], ydata=[height - handlebox.ydescent + LARGE_ADJUST_Y, 0 + LARGE_ADJUST_Y], linewidth=LARGE_SLASHWIDTH, color='black')
        handlebox.add_artist(l1)
        l2 = Line2D(xdata=[0 + LARGE_ADJUST_X, width - handlebox.xdescent], ydata=[0 + LARGE_ADJUST_Y, height - handlebox.ydescent + LARGE_ADJUST_Y], linewidth=LARGE_SLASHWIDTH, color='black')
        handlebox.add_artist(l2)
        return p
    
class LargeDotRectangleHandler(HandlerPatch):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        lower_left = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        center = 0.5 * (width - handlebox.xdescent + LARGE_ADJUST_X), 0.5 * (height - handlebox.ydescent + LARGE_ADJUST_Y)
        p = Rectangle(xy=lower_left, width=width, height=height, facecolor='white', edgecolor='black', linewidth=LARGE_LINEWIDTH)
        handlebox.add_artist(p)
        c = Circle(xy=center, radius=0.08 * width, facecolor='black')
        handlebox.add_artist(c)
        return p
    
custom_handler_map = {BackSlashRectangle: BackSlashHandler(), ForwardSlashRectangle: ForwardSlashHandler(), XRectangle: XHandler(), DotRectangle: DotRectangleHandler()}
large_custom_handler_map = {BackSlashRectangle: LargeBackSlashHandler(), ForwardSlashRectangle: LargeForwardSlashHandler(), XRectangle: LargeXHandler(), DotRectangle: LargeDotRectangleHandler()}