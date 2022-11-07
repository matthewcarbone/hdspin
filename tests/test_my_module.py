import hdspin
print("hdspin is", hdspin)

from hdspin.my_module import add


def test_add():
    assert add(1, 2) == 3
