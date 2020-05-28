import ipywidgets as wd

class function_box:
    def __init__(self, label, ctrl_dict, **kwargs):
        self.box_part = []
        self.slider_obj = {}
        # self.label = wd.Label(label)
        self.box_part.append(wd.Label(label))
        for var in kwargs:
            self.slider_obj.update({var : wd.IntSlider(value = kwargs[var][0], min = kwargs[var][1], max = kwargs[var][2], step = kwargs[var][3],
                            description= var + ':',
                            continuous_update=False)})
            self.box_part.append(self.slider_obj[var])
            ctrl_dict.update({var : self.slider_obj[var]})
        self.box = wd.VBox(children = self.box_part)

    def get_instance(self):
        return self.box
    def get_part(self):
        return self.box_part