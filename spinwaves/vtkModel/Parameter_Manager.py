class ParamManager():
    def __init__(self):
        self.parameters = []
        #each group of tied parameters will be given a group number
        self.groups = []#list of tuples (num(=index), bool(active/inactive))
    
    def AssignNewGroup(self, param):
        for group in self.groups:
            if not group[1]:
                param.group = group[0]
                group[1] = True
                return
        num = len(self.groups)
        self.groups.append( [num,True] )
        param.group = num
        #print "assigned ", param.group
    
    def addParam(self, param):
        """Appends a new parameter to this list."""
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
        param2.group = param1.group
        param2.default = param1.default
        #param1.fit = param2.fit
        #param1.value = param2.value
        #param1.min = param2.min
        #param1.max = param2.max
        #param1.group = param2.group
        
        self.__updateGroups()
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
    
    def __updateGroups(self):
        groups = []
        for param in self.parameters:
            for group in groups:
                if group == param.group:
                    break
            else:#Creating a list of unique groups that are represented
                groups.append(param.group)
        
        for i in range(len(self.groups)):
            for rep_group in groups:
                if rep_group == i:
                    break
            else:
                self.groups[i][1] = False
                
        #shuffle groups down to eliminate empty group numbers
        #for i in range(len(self.groups)):
        #    if not self.groups[i][1]:
        #        for param in self.parameters:
        #            if param.group > i:
        #                param.group -= 1
        #    self.groups[i][1] = True
                
    def getGroups(self):
        """Returns a list of lists of like-grouped parameters."""
        
    
    def untie(self, paramObj, index):
        """will untie the given Jparam object with the parameter given by index as well as
        all parameters tied to the parameter at index."""
        self.AssignNewGroup(paramObj)
        if paramObj.isTiedTo(index):
            paramObj.tied.remove(index)
            self.parameters[index].tied.remove(self.parameters.index(paramObj))
        
            for p in self.parameters[index].tied:
                self.untie(paramObj, p)
        self.__updateGroups()
    
    def removeParam(self, param):
        """removes the given JParam object from the list of parameters and corrects all
        indices and tied lists."""
        
        index = self.parameters.index(param)
        self.parameters.pop(index)
        for parameter in self.parameters:
            for tiedIndex in parameter.tied:
                if tiedIndex>=index:
                    tiedIndex-=1
        self.__updateGroups()
                    
    def getIndex(self, paramObj):
        """Returns the index of the object in the parameter list, or an error will occur if the
        parameter object is not in this manager."""
        return self.parameters.index(paramObj)
            