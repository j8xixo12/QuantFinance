import ipywidgets as wd

class function_box:
    def __init__(self, label,**kwargs):
        self.box_part = []
        # self.label = wd.Label(label)
        self.box_part.append(wd.Label(label))
        for var in kwargs:
            self.box_part.append(wd.IntSlider(value = kwargs[var][0], min = kwargs[var][1], max = kwargs[var][2], step = kwargs[var][3],
                            description= var + ':',
                            continuous_update=False))
        self.box = wd.VBox(children = self.box_part)

    def get_instance(self):
        return self.box
    def get_part(self):
        return self.box_part