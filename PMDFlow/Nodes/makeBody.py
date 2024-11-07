## Copyright 2024 - 2025 [Giacomo Cangi]

## Licensed under the MIT License;
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at

##                  https://mit-license.org/

## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.

#* standard packages
from PyFlow.Core import NodeBase
from PyFlow.Core.NodeBase import NodePinsSuggestionsHelper
from PyFlow import findPinClassByType, getAllPinClasses
from PyFlow import CreateRawPin
from PyFlow.Core.Common import *

#* my
from PMD.maker import Body

class makeBody(NodeBase):
    def __init__(self, name):
        super(makeBody, self).__init__(name)
        
        # Create Special pins
        # self.outexec = self.createOutputPin(pinName="outExec", dataType="ExecPin")
        # self.inexec = self.createInputPin(pinName="inExec", dataType="ExecPin", callback=self.processNode) 
        
        # Create Standard Input pins
        self.m = self.createInputPin("m (body mass) [kg]",'FloatPin', defaultValue=0.0, structure=StructureType.Multi, constraint="1", structConstraint="1")
        self.J = self.createInputPin("J (body inertia moment [kg*m^2])",'FloatPin', defaultValue=0.0, structure=StructureType.Multi, constraint="1", structConstraint="1")
        self.r = self.createInputPin("r (body CoG position vector [mm])",'FloatPin', structure=StructureType.Array, constraint="1", structConstraint="1")
        self.p = self.createInputPin("p (body orientation angle [rad])",'FloatPin', defaultValue=0.0, structure=StructureType.Multi, constraint="1", structConstraint="1")
        self.r_d = self.createInputPin("rd (body CoG position vector der., opt.) [m/s]",'FloatPin', structure=StructureType.Array, constraint="1", structConstraint="1")
        self.p_d = self.createInputPin("pd (body orientation vector der., opt) [rad/s]",'FloatPin', defaultValue=0.0, structure=StructureType.Multi, constraint="1", structConstraint="1")
        
        # Create Standard Output pins
        self.b = self.createOutputPin(pinName="body", dataType="InstancePin", defaultValue=None, structure=StructureType.Array, constraint=None, structConstraint=None, supportedPinDataTypes=[], group="")
        
    @staticmethod
    def pinTypeHints():
        helper = NodePinsSuggestionsHelper()
        helper.addInputDataType('AnyPin')
        helper.addOutputDataType('AnyPin')
        helper.addInputStruct(StructureType.Single)
        helper.addOutputStruct(StructureType.Single)
        return helper

    @staticmethod
    def category():
        return 'Multi-Body Maker'

    @staticmethod
    def keywords():
        return ["Make","Body"]

    @staticmethod
    def description():
        return "Create and add a body to the planar multi-body simulation."

    def compute(self, *args, **kwargs):
        
        #* Core - Body definition
        body = Body()
        body.m = self.m
        body.J = self.J 
        body.r = self.r (* 10**(-3)) # brings unit to m instad of mm
        body.p = self.p
        body.r_d = self.r_d
        body.p_d = self.p_d
        
        #* Send data to output
        self.b.setData(body)
        print(f"Body for simulation created!")