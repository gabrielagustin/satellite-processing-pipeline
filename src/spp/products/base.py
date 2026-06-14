from abc import ABC, abstractmethod

class ProductWriter(ABC):

    @abstractmethod
    def write(self, product):
        pass
