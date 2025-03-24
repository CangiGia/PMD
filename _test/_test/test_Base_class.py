import pytest
import numpy as np
from src.maker import Base

def test_instance_creation():
    # create instances of Base
    instance1 = Base()
    instance2 = Base()
    
    # check if instances are created
    assert isinstance(instance1, Base)
    assert isinstance(instance2, Base)

def test_instance_count():
    # reset the count for testing purposes
    Base.COUNT = 0
    
    # create instances of Base
    instance1 = Base()
    instance2 = Base()
    
    # check if the count is correct
    assert Base.get_count() == 2

def test_get_type():
    # create an instance of Base
    instance = Base()
    
    # check if the type is correct
    assert instance.get_type() == 'Base'

if __name__ == "__main__":
    pytest.main()