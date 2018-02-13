from random import normalvariate
from random import seed
from math import sqrt
#import matplotlib.pyplot as plt

class Node():
    numNodes = 0
    def __init__(self, x, y):
        self.isExitNode = False
        self.x = x
        self.y = y
        Node.numNodes += 1
        self.num = Node.numNodes
        self.riskA = {}
        self.riskB = {} 
        self.pointer = 0
        self.demand = 0

class Arc():
    numArcs = 0
    def __init__(self, head, tail, ffs, cap):
        Arc.numArcs += 1
        self.head = head
        self.tail = tail
        self.ffs = ffs     
        self.cap = cap
        self.fftt = 0
        distance = sqrt((self.head.x - self.tail.x)**2 + (self.head.y - self.tail.y)**2)
        if not self.ffs:
            self.fftt = 0
        elif self.ffs == 999999:
            self.fftt = 999999
        else:
            self.fftt = round(60 * distance / self.ffs, 8)

    def setFftt(self):
        distance = sqrt((self.head.x - self.tail.x)**2 + (self.head.y - self.tail.y)**2)
        if not self.ffs:
            self.fftt = 0
        else:
            self.fftt = round(60 * distance / self.ffs, 8)
        
class Grid(): 
    def __init__(self, gridSize):
        self.size = gridSize
        self.isPerturbed = False
        self.nodeArray = [[Node(x + 1, y + 1) for x in range(self.size)] for y in range(self.size)]
        self.arcList = []
        self.setExitNodes()

    def setExitNodes(self):
        for row in self.nodeArray:
            for node in row:
                if node.num < self.size + 1 or node.num % self.size == 1:
                    node.isExitNode = True

    def perturbGrid(self, mean, stdDev):
        if self.isPerturbed == False:
            for row in self.nodeArray:
                for node in row:
                    node.x += normalvariate(mean,stdDev)
                    node.y += normalvariate(mean,stdDev)
            self.isPerturbed = True

    def setArcs(self):
        n = len(self.nodeArray)
        for row in self.nodeArray:
            y = self.nodeArray.index(row)
            for node in row:
                x = row.index(node)
                if y > 0:
                    self.arcList.append(Arc(self.nodeArray[y - 1][x], node, 30, 200))
                if x > 0:
                    self.arcList.append(Arc(self.nodeArray[y][x - 1], node, 30, 200))
                if x < n - 1:
                    self.arcList.append(Arc(self.nodeArray[y][x + 1], node, 30, 200))
                if y < n - 1:
                    self.arcList.append(Arc(self.nodeArray[y + 1][x], node, 30, 200))

    def setRisks(self, scenario):
        for row in self.nodeArray:
            for node in row:
                node.riskA[scenario.num] = round(scenario.riskFuncA.getRisk(node.x, node.y), 2)
                node.riskB[scenario.num] = round(scenario.riskFuncB.getRisk(node.x, node.y), 2)

    def setDemand(self, mean, stdDev):
        for row in self.nodeArray:
            for node in row:
                if not node.isExitNode:
                    node.demand += round(normalvariate(mean,stdDev))
        
    def printGrid(self):
        for row in self.nodeArray:
            for node in row:
                print(node.num, node.demand, node.x, node.y, node.riskA, node.riskB)
        for arc in self.arcList:
            print(arc.head.num, arc.tail.num, arc.fftt, arc.cap)

    def plotGrid(self):
        xCoords = []
        yCoords = []
        for row in self.nodeArray:
            for node in row:
                xCoords.append(node.x)
                yCoords.append(-1 * node.y)
        plt.plot(xCoords, yCoords, 'ro') 
        plt.show()

class RiskFunction():
    def __init__(self, x, y, a, b, c, d, e, f):
        self.x = x
        self.y = y
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.f = f

    def getRisk(self, x, y):
        return (self.a*(x - self.x)**2
                + self.b*(y - self.y )**2
                + self.c*(x - self.x)*(y + self.y)
                + self.d*(x - self.x)
                + self.e*(y - self.y)
                + self.f
                )
                
class Scenario():
    _registry = []
    numScenarios = 0
    def __init__(self, probability, riskFuncA, riskFuncB):
        self.probability = probability
        self.riskFuncA = riskFuncA
        self.riskFuncB = riskFuncB
        Scenario.numScenarios += 1
        self.num = Scenario.numScenarios
        self._registry.append(self)

class Option():
    _registry = []
    numOptions = 0
    def __init__(self):
        Option.numOptions += 1
        self.num = Option.numOptions
        self._registry.append(self)
        self.accessArcList = [] 
        self.cost = 0
        self.arcCapacity = 0
        self.capacity = 0
        self.exclusionList = []
        self.isShelter = False

    def setShelterDetails(self, cost, capacity, node):
        self.isShelter = True
        self.cost = cost
        self.capacity = capacity
        self.node = node

    def setContraflowDetails(self, network, cost, ffs, capacity, entranceNodes, exitNodes, links):
        self.isShelter = False
        self.cost = cost
        self.arcCapacity = capacity
        self.entranceNodes = entranceNodes
        self.exitNodes = exitNodes
        self.allNodes = list(set(entranceNodes) | set(exitNodes))
        #print([node.num for node in self.allNodes])

        #duplicate and connect nodes
        for node in self.allNodes:
            newNode = Node(node.x, node.y)
            for link in links:
                for accessNode in link:
                    if accessNode.num == node.num:
                        link[link.index(accessNode)] = newNode 
            network.nodeList.append(newNode)
            if node in entranceNodes:
                network.arcList.append(Arc(newNode, node, 0, 1))
                self.accessArcList.append(network.arcList[-1])
            else:
                network.arcList.append(Arc(newNode, node, 999999, 1))
                self.accessArcList.append(network.arcList[-1])
            if node in exitNodes:
                network.arcList.append(Arc(node, newNode, 0, 1))
                self.accessArcList.append(network.arcList[-1])
            else:
                network.arcList.append(Arc(node, newNode, 999999, 1))
                self.accessArcList.append(network.arcList[-1])

        #add real arcs
        for link in links:
            network.addRealArcPair(link[0], link[1], ffs, capacity)

class RedpNetwork():
    def __init__(self):
        self.name = '' 
        self.inputFilePointer = []
        self.comment = ''
        self.grid = []
        self.nodeList = []
        self.arcList = []
        self.budget = 0
        self.sinkNode = []

    def setBudget(self, budget):
        self.budget = budget

    def setSinkNode(self):
        self.sinkNode = Node.numNodes + 1

    def sortArcList(self):
        self.arcList.sort(key = lambda x: x.tail.num)

    def setNodePointers(self):
        for node in self.nodeList:
            node.pointer = Arc.numArcs
        for arc in self.arcList:
            if self.arcList.index(arc) + 1 < self.nodeList[arc.tail.num - 1].pointer:
                self.nodeList[arc.tail.num - 1].pointer = self.arcList.index(arc) + 1

    def setInputFilePointer(self):
        self.inputFilePointer.append(3)
        self.inputFilePointer.append(self.inputFilePointer[-1] + Node.numNodes)
        self.inputFilePointer.append(self.inputFilePointer[-1] + Arc.numArcs)
        self.inputFilePointer.append(self.inputFilePointer[-1] + Option.numOptions)
        self.inputFilePointer.append(self.inputFilePointer[-1] + Option.numOptions)
        self.inputFilePointer.append(self.inputFilePointer[-1] + Option.numOptions)
        self.inputFilePointer.append(self.inputFilePointer[-1] + Option.numOptions)
        self.inputFilePointer.append(self.inputFilePointer[-1] + Scenario.numScenarios)
        self.inputFilePointer.append(self.inputFilePointer[-1] + Scenario.numScenarios)
        self.inputFilePointer.append(self.inputFilePointer[-1] + Scenario.numScenarios)
        self.inputFilePointer.append(self.inputFilePointer[-1] + 1)

    def addAbstractArcPair(self, node1, node2, isOpenTo1, isOpenTo2):
        self.arcList.append(Arc(node1, node2, 0 if isOpenTo1 else 999999, 1))
        self.arcList.append(Arc(node2, node1, 0 if isOpenTo2 else 999999, 1))

    def addRealArcPair(self, node1, node2, ffs, cap):
        self.arcList.append(Arc(node1, node2, ffs, cap))
        self.arcList.append(Arc(node2, node1, ffs, cap))

    def constructOptions(self):
        #create do-nothing case
        self.doNothingOption = Option()
        innerNode = self.nodeList[0]
        middleNode = self.nodeList[(self.grid.size**2-1)//2]
        outerNode = self.nodeList[self.grid.size**2-1]
        self.doNothingOption.setContraflowDetails(self, 0, 60, 200,
            [innerNode, outerNode, middleNode],
            [innerNode, outerNode, middleNode],
            [[innerNode, middleNode], [middleNode, outerNode]])

        #create contraflow case
        self.contraflowOption = Option()
        self.contraflowOption.setContraflowDetails(self, 8000, 60, 400,
            [outerNode, middleNode],
            [innerNode],
            [[innerNode, middleNode], [middleNode, outerNode]])

        #make options mutually exclusive
        self.doNothingOption.exclusionList.append(self.contraflowOption)
        self.contraflowOption.exclusionList.append(self.doNothingOption)

        #add shelter options
        self.shelterOption1 = Option()
        self.shelterOption1.setShelterDetails(4000, 5000, self.nodeList[(self.grid.size**2-1)//3])
        self.shelterOption2 = Option()
        self.shelterOption2.setShelterDetails(4000, 5000, self.nodeList[2*(self.grid.size**2-1)//3])
        self.shelterOption3 = Option()
        self.shelterOption3.setShelterDetails(4000, 5000, self.nodeList[3*(self.grid.size**2-1)//4])

        #add sink node
        self.sinkNode = Node(0,0)
        self.nodeList.append(self.sinkNode)

        #add shelter and shelter in place arcs
        for row in self.grid.nodeArray:
            for node in row:
                newArc = Arc(self.sinkNode, node, 0, 1)
                self.arcList.append(newArc)
                for option in Option._registry:
                    if option.isShelter and node == option.node:
                        option.accessArcList = [newArc]
                self.arcList.append(Arc(node, self.sinkNode, 0, 1))

    def constructNetwork(self):
        self.constructOptions()
        self.sortArcList()
        self.setNodePointers()
        self.setInputFilePointer()
        self.gridSizeName = str(self.grid.size) + 'X' + str(self.grid.size) 
        self.name = self.gridSizeName + 'Test'
        self.comment = 'Network to test generation script with ' + self.gridSizeName + ' grid, ' + str(Scenario.numScenarios) + ' scenarios and ' + str(Option.numOptions) + ' options.' 

    def setGrid(self, grid):
        self.grid = grid
        for row in self.grid.nodeArray:
            for node in row:
                self.nodeList.append(node)
        for arc in self.grid.arcList:
            self.arcList.append(arc)

    def addOption(self, option):
        self.options.append[option]

    def printNetwork(self):
        print(self.name, " ".join(str(pointer) for pointer in self.inputFilePointer))
        print(self.comment)
        for node in self.nodeList:
            print(node.pointer, node.demand)
        for arc in self.arcList:
            print(arc.head.num, arc.tail.num, arc.fftt, arc.cap)
        for option in Option._registry:
            print(" ".join(str(self.arcList.index(arc)) for arc in option.accessArcList))
        for option in Option._registry:
            if option.exclusionList:
                print(" ".join(str(exclusion.num) for exclusion in option.exclusionList), end = ' ')
                print()
            else:
                print('0')   
        for option in Option._registry:
            print(option.capacity)
        for option in Option._registry:
            print(option.cost)
        for scenario in Scenario._registry:
            print(scenario.probability)
        for i in range(Scenario.numScenarios):
            for row in self.grid.nodeArray:
                for node in row:
                    print(node.riskA[i + 1], end=' ')
            print()
        for i in range(Scenario.numScenarios):
            for row in self.grid.nodeArray:
                for node in row:
                    print(node.riskB[i + 1], end=' ')
            print()
        print(self.budget)
        print(self.sinkNode.num)

    def printNetworkFile(self):
        f = open(self.name + '.txt', 'w')
        print(self.name, " ".join(str(pointer) for pointer in self.inputFilePointer), file = f)
        print(self.comment, file = f)
        for node in self.nodeList:
            print(node.pointer, node.demand, file = f)
        for arc in self.arcList:
            print(arc.head.num, arc.tail.num, arc.fftt, arc.cap, file = f)
        for option in Option._registry:
            print(" ".join(str(self.arcList.index(arc)) for arc in option.accessArcList), end = ' ', file = f)
            print(file = f)
        for option in Option._registry:
            if option.exclusionList:
                print(" ".join(str(exclusion.num) for exclusion in option.exclusionList), end = ' ', file=f)
                print(file = f)
            else:
                print('0', file = f)   
        for option in Option._registry:
            print(option.capacity, file = f)
        for option in Option._registry:
            print(option.cost, file=f)
        for scenario in Scenario._registry:
            print(scenario.probability, file=f)
        for i in range(Scenario.numScenarios):
            for row in self.grid.nodeArray:
                for node in row:
                    print(node.riskA[i + 1], end=' ', file=f)
            print(file = f)
        for i in range(Scenario.numScenarios):
            for row in self.grid.nodeArray:
                for node in row:
                    print(node.riskB[i + 1], end=' ', file=f)
            print(file = f)
        print(self.budget, file=f)
        print(self.sinkNode.num, file=f)
        
seed(0)
n = 3
gridSize = 2**n + 1
gridMean = 0
gridStdDev = 0.1
demandMean = 500
demandStdDev = 50
budget = 12000
maxRisk = 10
riskSpread = -maxRisk/2/(gridSize + 1)**2

#test
testRisk1 = RiskFunction(gridSize + 1, gridSize + 1, riskSpread, riskSpread, 0, 0, 0, maxRisk)
testRisk2 = RiskFunction(0, 0, riskSpread, riskSpread, 0, 0, 0, maxRisk)
testRisk3 = RiskFunction(0, gridSize + 1, riskSpread, riskSpread, 0, 0, 0, maxRisk)
testRisk4 = RiskFunction(gridSize + 1, 0, riskSpread, riskSpread, 0, 0, 0, maxRisk)
testScenario1 = Scenario(0.10, testRisk1, testRisk2)
testScenario2 = Scenario(0.20, testRisk2, testRisk2)
testScenario3 = Scenario(0.30, testRisk3, testRisk2)
testScenario4 = Scenario(0.40, testRisk4, testRisk2)
testGrid = Grid(gridSize)
testGrid.perturbGrid(gridMean, gridStdDev)
testGrid.setRisks(testScenario1)
testGrid.setRisks(testScenario2)
testGrid.setRisks(testScenario3)
testGrid.setRisks(testScenario4)
testGrid.setDemand(demandMean, demandStdDev)
testGrid.setArcs()
#testOption1 = ContraflowOption(testGrid.nodeArray[0][0],
        #testGrid.nodeArray[2][2], 4000, 60, 200)
#testOption2 = ContraflowOption(testGrid.nodeArray[0][0],
        #testGrid.nodeArray[2][2], 4000, 60, 200)
#testGrid.plotGrid()
testNetwork = RedpNetwork()
testNetwork.setGrid(testGrid)
testNetwork.setBudget(budget)
testNetwork.constructNetwork()
testNetwork.printNetwork()
testNetwork.printNetworkFile()
