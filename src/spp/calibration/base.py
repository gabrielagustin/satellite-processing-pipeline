from abc import ABC, abstractmethod

class Calibrator(ABC):

    @abstractmethod
    def apply(self, image):
        pass
