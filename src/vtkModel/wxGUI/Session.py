import xml.dom.minidom
import xml.dom.ext

class Session():
    """Stores information about a user session"""
    def __init__(self):
        self.bondTable = None
        self.atomTable = None
        
    def saveToXML(self, filename):
        #Create the document
        doc = xml.dom.minidom.Document()
        
        #Create the main element
        main_element = doc.createElement('Spinwaves Session')
        doc.appendChild(main_element)
        
        element2 = doc.createElement('element')
        main_element.appendChild(element2)
        element2.setAttribute('name', 'value')
        description = doc.createTextNode("A quiet, scenic park with lots of wildlife.")
        element2.appendChild(description)
        
        
        
        #Write to screen
        xml.dom.ext.PrettyPrint(doc)
        
        #Write to the file
        xml.dom.ext.PrettyPrint(doc, open(filename, 'w'))