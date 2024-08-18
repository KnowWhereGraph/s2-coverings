from ..src.lib.config import Config


def test_defaults():
    """
    The default values were carefully chosen - test to make sure on one has changed them

    :return: None
    """
    config = Config()
    assert config.tolerance == 1e-2
    assert config.min_level == 13
    assert config.max_level == 13
