# Generated by Django 3.1.4 on 2021-03-12 09:26

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('wooey', '0043_update_model_protection'),
    ]

    operations = [
        migrations.AlterField(
            model_name='scriptversion',
            name='script_path',
            field=models.FileField(upload_to=''),
        ),
        migrations.AlterField(
            model_name='wooeyfile',
            name='filepath',
            field=models.FileField(max_length=500, upload_to=''),
        ),
    ]
