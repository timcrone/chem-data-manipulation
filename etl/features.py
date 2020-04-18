"""
This library contains tools manipulating features in the
ZINC database.

Features are one or several of:
  aggregator
  biogenic
  endogenous
  fda
  in-man
  investigational
  in-vitro
  in-vivo
  metabolite
  world

This is stored as a bitmap for convenience.
"""

import logging


class Features:
    """
    Defines methods for working with features and feature maps.

    The class is lazy, pass in an initialization list or map
    and the corresponding string or map representation will
    be calculated at the necessary time.
    """
    KNOWN_FEATURES = {
        'aggregator': 1,
        'biogenic': 2,
        'endogenous': 4,
        'fda': 8,
        'in-man': 16,
        'investigational': 32,
        'in-vitro': 64,
        'in-vivo': 128,
        'metabolite': 256,
        'world': 512
    }

    def __init__(self, feature_str=None, feature_map=None):
        """
        Pass the comma-delimited feature string
        or the feature bitmap to initialize.
        feature_list takes precedence if passed.

        Feature string is available in .str()
        Feature map is available in .map()

        :param feature_str: list of features in comma-delimited form
        :type feature_str: str
        :param feature_map: bitmap of features
        :type feature_map: int
        """
        self.feature_map = None
        self.feature_str = None
        if feature_str is not None:
            self.feature_str = feature_str
        elif feature_map is not None:
            self.feature_map = feature_map

    def map(self):
        """
        Returns the mapped value of the feature list.

        :return: bitmap of features in KNOWN_FEATURES
        :rtype: int
        """
        if self.feature_map is not None:
            pass
        elif self.feature_str is not None:
            self.feature_map = self._features(feature_list=self.feature_str)
        else:
            raise ValueError("Initialize with feature string or map")
        return self.feature_map

    def str(self):
        """
        Returns the string value of the feature list.

        :return: comma-separated string of features
        :rtype: str
        """
        if self.feature_str is not None:
            pass
        elif self.feature_map is not None:
            self.feature_str = self._features(feature_map=self.feature_map)
        else:
            raise ValueError("Initialize with feature string or map")
        return self.feature_str

    def _features(self, feature_list=None, feature_map=None):
        """
        Map features to bits or vice versa.  If feature_list is passed
        then the bitmap is returned, if feature_map is passed then the
        string form is returned.
        :param feature_list: list of features in comma-delimited form
        :type feature_list: str
        :param feature_map: bitmap of features
        :type feature_map: int
        :return: bitmap or list of features
        :rtype: int|str
        """
        if feature_list is not None:
            _result = 0
            for feature in feature_list.split(','):
                if feature:
                    _result |= self.KNOWN_FEATURES[feature]
        elif feature_map is not None:
            _result = ""
            for feature in self.KNOWN_FEATURES:
                if feature and feature_map & self.KNOWN_FEATURES[feature]:
                    _result += f"{feature},"
            _result = _result[:-1]
        else:
            return 0
        return _result


if __name__ == "__main__":
    # pylint: disable=invalid-name
    logging.basicConfig(level=logging.INFO)
    test_string = "in-vivo,in-vitro"
    test = Features(feature_str=test_string)
    assert test.map() == 64 | 128
    assert test.str() == test_string
    test_map = 512 | 16 | 8
    test = Features(feature_map=test_map)
    assert test.map() == test_map
    assert test.str() == "fda,in-man,world"
    test_map = 1023
    test = Features(feature_map=test_map)
    assert test.str().split(',') == list(test.KNOWN_FEATURES.keys())
    test = Features(feature_str=','.join((test.KNOWN_FEATURES.keys())))
    assert test.map() == 1023
    test = Features(feature_map=0)
    assert test.str() == ''
    assert test.map() == 0
    test = Features(feature_str='')
    assert test.map() == 0
    assert test.str() == ''
