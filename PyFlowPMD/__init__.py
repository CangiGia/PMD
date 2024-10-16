PACKAGE_NAME = 'PyFlowPMD'
import os
from collections import OrderedDict
from PyFlow.UI.UIInterfaces import IPackage


#* Retrieve the directory containing this __init__.py file
master_folder_path = os.path.dirname(__file__)

#* Class based nodes - AUTO IMPORT FROM FOLDER
print(f"\n\nImporting Nodes from {PACKAGE_NAME}...")
node_folder_path = os.path.join(master_folder_path, 'Nodes')
_NODES = {}

# Get a list of all files in the folder
files = os.listdir(node_folder_path)

# Iterate over the files and import the modules
for file in files:
    if not file.startswith("_") and file.endswith(".py"): # Only .py files, files starting with "_" are skipped.
        module_name = file[:-3]  # Remove the ".py" extension
        print(f"Loading Node:\t{module_name}")
        
        import_statement = f"from PyFlow.Packages.PyFlowPMD.Nodes.{module_name} import {module_name}"
        exec(import_statement)

		# Load the nodes...
        _NODE_i = f"_NODES[{module_name}.__name__] = {module_name}"
        exec(_NODE_i)

#* Factories
_FOO_LIBS = {}
# _NODES = {} # no more required, replaced by the above routine
_PINS = {}
_TOOLS = OrderedDict()
_PREFS_WIDGETS = OrderedDict()
_EXPORTERS = OrderedDict()

class PyFlowPMD(IPackage):
	def __init__(self):
		super(PyFlowPMD, self).__init__()

	@staticmethod
	def GetExporters():
		return _EXPORTERS

	@staticmethod
	def GetFunctionLibraries():
		return _FOO_LIBS

	@staticmethod
	def GetNodeClasses():
		return _NODES

	@staticmethod
	def GetPinClasses():
		return _PINS

	@staticmethod
	def GetToolClasses():
		return _TOOLS