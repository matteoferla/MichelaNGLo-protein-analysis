## _BaseMixin contains __getattr__ and _failsafe decorator
from warnings import warn
from michelanglo_protein.settings_handler import global_settings #the instance not the class.


class _BaseMixin:
    settings = global_settings
    # these older class attributes should be made redundant
    #    fetch = True
    #    croak = True
    #    tollerate_no_SNV = True

    # decorator
    def _failsafe(func):
        def wrapper(self, *args, **kargs):
            # the call happned after chekcing if it should croak on error so to make the traceback cleaner.
            if self.settings.error_tolerant:
                try:
                    return func(self, *args, **kargs)
                except Exception as error:
                    print('Error caught in method `Protein().{n}`: {e}'.format(n=func.__name__, e=error))
                    return None
            else:
                return func(self, *args, **kargs)

        return wrapper


    def __getattr__(self, item):
        if self.settings.missing_attribute_tolerant:
            if item not in self.__dict__:
                if item in self.other:  ## it is in the trash!
                    warn('Accessed attribute in other list. Thanks for proving the key:value pair. But please dont abuse this backdoor!')
                    self.__dict__[item] = self.other[item]
                else:
                    warn('Accessed non-existant attribute {item} for Protein instance. Likely cause the code changed, but the from_pickle flag is True.'.format(v=self, item=item))
                    self.__dict__[item] = 'Unknown'
        return self.__dict__[item]

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, state):
        self.__dict__ = state
