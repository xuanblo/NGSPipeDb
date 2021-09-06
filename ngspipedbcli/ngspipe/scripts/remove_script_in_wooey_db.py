import sys
import os


workdir = sys.argv[1]

sys.path.append(workdir)
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'ngsdb.settings')

import django
django.setup()

from wooey.models.core import Script

scripts = Script.objects.all()
print('before:')
print(scripts)
print(len(scripts))

scripts = Script.objects.all().delete()

scripts = Script.objects.all()
print('after delete')
print(scripts)
print(len(scripts))
