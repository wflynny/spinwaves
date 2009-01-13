class ParamManager():
    def __init__(self):
        self.parameters = []
    
    def addParam(self, param):
        """Adds a new parameter to this list."""
        print "addParam function"
        self.parameters.append(param)
        
    def validIndex(self, index):
        """Checks if the integer index is a valid index(identifier) of a parameter."""
        return (index < len(self.parameters))
    
    def tie(self, paramObj, param2_index):
        """Adds param2_index to the list of tied indices on paramObj, and vice versa.
        Also changes the min, max, value, and fit(bool) values in the parameter given
        by param2_index to equal those of paramObj."""
        index1 = self.parameters.index(paramObj)
        param1 = paramObj
        index2 = param2_index
        param2 = self.parameters[index2]
        if index1 != index2:
            #Add second index to first param
            if not param1.isTiedTo(index2):
                param1.tied.append(index2)
            #Add first index to second param
            if not param2.isTiedTo(index1):
                param2.tied.append(index1)
            
        #Change values on second parameter to equal those of paramObj
        param2.fit = param1.fit
        param2.value = param1.value
        param2.min = param1.min
        param2.max = param1.max
        
        #tie to any parameters tied to the second parameter
        for eachIndex in param2.tied:
            if not param1.isTiedTo(eachIndex):
                self.tie(param1, eachIndex)
        
        for eachIndex in param1.tied:
            if not param2.isTiedTo(eachIndex):
                self.tie(param2, eachIndex)
                
        #to look nicer for the user
        param1.tied.sort()
        param2.tied.sort()
    
    def untie(self, paramObj, index):
        paramObj.tied.remove(index)
        self.parameters[index].tied.remove(self.parameters.index(paramObj))
    
    def removeParam(self, param):
        """removes the given JParam object from the list of parameters and corrects all
        indices and tied lists."""
        
        index = self.parameters.index(param)
        self.parameters.pop(index)
        for parameter in self.parameters:
            for tiedIndex in parameter.tied:
                if tiedIndex>=index:
                    tiedIndex-=1
                    
    def getIndex(self, paramObj):
        """Returns the index of the object in the parameter list, or an error will occur if the
        parameter object is not in this manager."""
        return self.parameters.index(paramObj)
            