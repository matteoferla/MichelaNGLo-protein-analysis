import os

def get_default_params_filenames():
    """
    This was never meant to stay for long...
    I was hoping to switch to having the params passed as strings as oppossed to written to file!
    When that is corrected, this will use `importlib.resources.read_text`

    :return:
    """
    default_params_folder = os.path.dirname(__file__)
    return [os.path.join(default_params_folder, filename) for filename in os.listdir(default_params_folder)
                      if os.path.splitext(filename)[1] == '.params']

