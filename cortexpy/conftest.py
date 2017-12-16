from hypothesis import settings

settings.register_profile("ci", settings(max_examples=1000))
settings.register_profile("dev", settings(max_examples=5))
# settings.register_profile("debug", settings(max_examples=10, verbosity=Verbosity.verbose))
# settings.load_profile(os.getenv('HYPOTHESIS_PROFILE', 'default'))
