from abc import ABC, abstractmethod


class QAValidator(ABC):
    """
    Abstract interface for product quality validation.

    Implementations are responsible for evaluating the quality
    and consistency of generated products.

    Typical checks may include:

    - Invalid pixel detection
    - Saturation detection
    - Metadata consistency
    - Radiometric sanity checks
    - Product completeness verification

    Validators should not modify products. Their role is limited
    to assessment and reporting.
    """

    @abstractmethod
    def validate(self, product):
        """
        Validate a processing product.

        Parameters
        ----------
        product
            Product to be evaluated.

        Returns
        -------
        dict
            Validation results and quality metrics.
        """
        pass
