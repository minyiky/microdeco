from microdeco.engine import Engine
from microdeco.model import ZH_L16B_GF

model = ZH_L16B_GF(gf_low=0.30, gf_high=0.80)
engine = Engine(model=model)
engine.add_gas(21, depth=0)

profile = engine.calculate(35, 40)

for s in profile:
    print('Step(phase="{}", abs_p={:.4f}, time={:.4f},' \
    ' gf={:.4f})'.format(s.phase, s.abs_p, s.time, s.data.gf))

for stop in engine.deco_table:
    print(stop)

print(engine.deco_table.total)
