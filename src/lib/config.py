import os


class Config:
    max_level = int(os.getenv("MAX_LEVEL", 13))
    min_level = int(os.getenv("MIN_LEVEL", 13))
    tolerance = float(os.getenv("TOLERANCE", 1e-2))


config = Config()
