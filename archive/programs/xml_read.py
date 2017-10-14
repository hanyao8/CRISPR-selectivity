from xml.dom import minidom
Test_file = open('C:\Users\Choon\Downloads\sequence.xml','r')
xmldoc = minidom.parse(Test_file)

Test_file.close()

def printNode(node):
  print node
  for child in node.childNodes:
       printNode(child)

printNode(xmldoc.documentElement)