class Component:
    def __init__(self, name):
        self.name = name
        self.inputs = []
        self.outputs = []

    def add_input(self, connection):
        self.inputs.append(connection)

    def add_output(self, connection):
        self.outputs.append(connection)

class CycleCloser(Component):
    def __init__(self, name):
        super().__init__(name)

class Compressor(Component):
    def __init__(self, name):
        super().__init__(name)

class Valve(Component):
    def __init__(self, name):
        super().__init__(name)

class SimpleHeatExchanger(Component):
    def __init__(self, name, heat_exchanger_type):
        super().__init__(name)
        self.heat_exchanger_type = heat_exchanger_type

# Define Connection class separately
class Connection:
    def __init__(self, source_component, source_port, target_component, target_port, label):
        self.source_component = source_component
        self.source_port = source_port
        self.target_component = target_component
        self.target_port = target_port
        self.label = label

# Example usage:
cc = CycleCloser('cycle closer')
cp = Compressor('compressor')
va = Valve('expansion valve')
co = SimpleHeatExchanger('condenser', 'condenser')
ev = SimpleHeatExchanger('evaporator', 'evaporator')

c1 = Connection(cc, 'out1', ev, 'in1', label='1')
c2 = Connection(ev, 'out1', cp, 'in1', label='2')
c3 = Connection(cp, 'out1', co, 'in1', label='3')
c4 = Connection(co, 'out1', va, 'in1', label='4')
c0 = Connection(va, 'out1', cc, 'in1', label='0')

co.add_input(c3)
co.add_output(c4)
ev.add_input(c1)
ev.add_output(c2)
cp.add_input(c2)
cp.add_output(c3)
va.add_input(c4)
va.add_output(c0)

# Accessing inputs and outputs of a component
print("Inputs of Simple Heat Exchanger 'condenser':", [conn.label for conn in co.inputs])
print("Outputs of Simple Heat Exchanger 'condenser':", [conn.label for conn in co.outputs])
