class Simulation(object):
	def __init__(self, size=2e-3, resolution = 50,model = None, simType='2d', mode='dots', debugger = True):
		super(Simulation, self).__init__()
		print("[!] Debugger is on by default. To turn it off specify \"debugger = False\" in the Simualtion object declaration"*debugger)
		self.debugger = debugger
		self.model, self.type, self.mode = model, simType, mode 
		self.start = lambda: self.model.plot(resolution, size)

	def setModel(self, newModel:str):
		if not newModel: raise TypeError("Can't set simulation model to None")
		if self.model == None: 
			self.model = newModel
			self.model.debugger = self.debugger
			return self
		raise ValueError(f"Simulation model already set to {self.model}")