from abc import ABC, abstractmethod


class ProcessingPipeline(ABC):
    """
    Abstract interface for processing pipelines.

    Pipelines orchestrate the execution of the different
    processing stages required to generate a satellite product.

    Typical stages include:

    - Data ingestion
    - Calibration
    - Quality assessment
    - Product generation

    The pipeline coordinates execution but should not contain
    sensor-specific processing logic.
    """

    @abstractmethod
    def run(self):
        """
        Execute the processing workflow.

        Returns
        -------
        object
            Generated processing product.
        """
        pass
