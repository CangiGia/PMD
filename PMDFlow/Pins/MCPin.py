from PyFlow.Core import PinBase
from PyFlow.Core.Common import *
import json

'''

Replace this FakeTypeRSPP with some custom class that can be created as FakeTypeRSSPP(value)

It can be used with numpy for example as this:

import numpy as np

@staticmethod
def pinDataTypeHint():
    return 'NumpyArray', np.zeros((2, 2, 3), np.uint8)

@staticmethod
def internalDataStructure():
    return np.ndarray
    
@staticmethod
def processData(data):
    if data is None:
        return NumpyArray.pinDataTypeHint()[1]
    if isinstance(data, np.ndarray):
        return data
    else:
        raise Exception("non Valid Image//numpy Array")
'''

class FakeTypeRSSPP(object):
    """docstring for FakeTypeRSSPP"""
    def __init__(self, value=None):
        super(FakeTypeRSSPP, self).__init__()
        self.value = value

class NoneEncoder(json.JSONEncoder):
    def default(self, vec3):
        return None

class NoneDecoder(json.JSONDecoder):
    def init(self, args, **kwargs):
        super(NoneDecoder, self).init(object_hook=self.object_hook, *args, **kwargs)

    def object_hook(self, vec3Dict):
        return None

class MCPin(PinBase):
    """doc string for MCPin"""
    def __init__(self, name, parent, direction, **kwargs):
        super(MCPin, self).__init__(name, parent, direction, **kwargs)
        self.setDefaultValue(MCPin.pinDataTypeHint()[1])
        self.disableOptions(PinOptions.Storable) #* DISABLE STORE CAPABILITY
    
    @staticmethod
    def jsonEncoderClass():
        return NoneEncoder

    @staticmethod
    def jsonDecoderClass():
        return NoneDecoder
    
    @staticmethod
    def IsValuePin():
        return True

    @staticmethod
    def supportedDataTypes():
        return ('MCPin',) 

    @staticmethod
    def pinDataTypeHint(): # This is datatype name and default Value
        return 'MCPin', FakeTypeRSSPP()

    @staticmethod
    def color():
        return (255, 165, 0, 255)

    @staticmethod
    def internalDataStructure():
        return FakeTypeRSSPP

    @staticmethod
    def processData(data):
        if data is None:
            return MCPin.pinDataTypeHint()[1]
        if isinstance(data, FakeTypeRSSPP):
            return data
        else:
            raise Exception("non Valid FakeTypeRSSPP")