from microdeco.models import ZH_L16C_GF, GradientFactors
from microdeco.gas import Gas


def ceil_to_stop(ceiling: float) -> int:
    """Round ceiling depth up to nearest 3 meters for standard deco stops."""
    if ceiling == 0:
        return 0
    return (int(ceiling // 3) + 1) * 3


if __name__ == "__main__":
    gas = Gas(p_o2=21, p_n2=79, p_he=0)
    gf = GradientFactors(gf_low=0.3, gf_high=0.85)
    model = ZH_L16C_GF()

    time = 0.0
    current_depth = 0
    step_duration = 0.1  # minutes

    # Phase 1: Descent to 30m at ~10 m/min (1 bar/min)
    print("=== DESCENT ===")
    descent_rate = 10  # meters/min
    target_depth = 30
    descent_time = target_depth / descent_rate
    num_steps = int(descent_time / step_duration)

    for i in range(num_steps):
        time += step_duration
        next_depth = (target_depth / descent_time) * time
        model.load_tissues_with_rate(current_depth, next_depth, step_duration, gas)
        ceiling = model.calculate_ceilings(gf)
        print(f"{time:6.1f} min (depth {next_depth:5.1f}m): ceiling = {ceiling:6.1f}m")
        current_depth = next_depth

    # Phase 2: Bottom time at 30m for 20 minutes
    print("\n=== BOTTOM TIME ===")
    bottom_duration = 20
    num_steps = int(bottom_duration / step_duration)

    for i in range(num_steps):
        time += step_duration
        model.load_tissues(target_depth, step_duration, gas)
        ceiling = model.calculate_ceilings(gf)
        print(f"{time:6.1f} min (depth {target_depth:5.1f}m): ceiling = {ceiling:6.1f}m")

    # Phase 3: Ascent from 30m to surface at ~10 m/min
    print("\n=== ASCENT ===")
    ascent_rate = 10  # meters/min
    ascent_distance = target_depth
    ascent_time = ascent_distance / ascent_rate
    num_steps = int(ascent_time / step_duration)

    for i in range(num_steps):
        time += step_duration
        next_depth = target_depth - (ascent_rate / ascent_time) * (time - (3 + bottom_duration))
        next_depth = max(next_depth, 0)
        model.load_tissues_with_rate(current_depth, next_depth, step_duration, gas)
        ceiling = model.calculate_ceilings(gf)
        stop_depth = ceil_to_stop(ceiling)
        print(f"{time:6.1f} min (depth {next_depth:5.1f}m): ceiling = {ceiling:6.1f}m, stop = {stop_depth:3.0f}m")
        current_depth = next_depth

    print("\n=== SURFACE ===")
    print(f"Total dive time: {time:.1f} minutes")
    print(f"Final ceiling: {model.calculate_ceilings(gf):.1f}m")
