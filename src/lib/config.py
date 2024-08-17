import os


class Config:
    max_level = os.getenv("MAXLEVEL", 13)
    min_level = os.getenv("MIN_LEVEL", 13)
    tolerance = os.getenv("TOLERANCE", 1e-2)


config = Config()
