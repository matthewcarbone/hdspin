from hdspin.simulation_engine import run_simulation


def test_sim_engine():
    sim = run_simulation(5, 6)
    assert sim.test_return == 5
    assert sim.test_return2 == 6
