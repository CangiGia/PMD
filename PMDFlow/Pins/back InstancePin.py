from PyFlow.Core import PinBase
from PyFlow.Core.Common import *
import json

'''

This Pin pass an object (class) instance.
pin shall be non-storable since it is not serializable.
'''
class FakeTypeRSSPP(object):        #? what is this doing?
    """docstring for FakeTypeRSSPP"""
    def __init__(self, value=None):
        super(FakeTypeRSSPP, self).__init__()
        self.value = value
        print("xxxxxxxxxxxxxxxxxxxxxxxxxxx") #! temp

class NoneEncoder(json.JSONEncoder): #!!!!!!
    def default(self, vec3):
        return None

class NoneDecoder(json.JSONDecoder): #!!!!!!
    def init(self, args, **kwargs):
        super(NoneDecoder, self).init(object_hook=self.object_hook, *args, **kwargs)

    def object_hook(self, vec3Dict):
        print("xxx Hook called...") #wip remove it
        return None #self.trimesh.base.Trimesh()    # I don't have trimesh #? What should I put here? 

class InstancePin(PinBase):
    """doc string for DemoPin"""
    def __init__(self, name, parent, direction, **kwargs):
        super(InstancePin, self).__init__(name, parent, direction, **kwargs)
        self.setDefaultValue(False)
        self.disableOptions(PinOptions.Storable) #* DISABLE STORE CAPABILITY

    @staticmethod
    def IsValuePin():
        return True

    # @staticmethod
    # def supportedDataTypes():
    #     return 'InstancePin, List' #! non so sicuro

    @staticmethod
    def pinDataTypeHint():
        return 'Instance of an object', False

    @staticmethod
    def color():
        return (200, 200, 50, 255)

    @staticmethod                       #? what is this doing?
    def internalDataStructure():
        return FakeTypeRSSPP
    
    @staticmethod
    def internalDataStructure():
        print("called internalDataStructure")
        return list #! what shoud I set? tried with object, float, list

    @staticmethod
    def processData(data):
        return InstancePin.internalDataStructure()(data)
    
    
    
    #! test from here
    @staticmethod
    def jsonEncoderClass():
        return NoneEncoder

    @staticmethod
    def jsonDecoderClass():
        return NoneDecoder