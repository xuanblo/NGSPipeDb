# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#   * Rearrange models' order
#   * Make sure each model has one field with primary_key=True
#   * Make sure each ForeignKey and OneToOneField has `on_delete` set to the desired behavior
#   * Remove `managed = False` lines if you wish to allow Django to create, modify, and delete the table
# Feel free to rename the models, but don't rename db_table values or field names.
from django.db import models


class Autoincrements(models.Model):
    base = models.TextField(blank=True, null=True)
    n = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'autoincrements'


class Directives(models.Model):
    directive = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'directives'


class Duplicates(models.Model):
    idspecid = models.TextField(blank=True, null=True)
    newid = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'duplicates'


class Features(models.Model):
    id = models.TextField(blank=True, null=False, primary_key=True)
    seqid = models.TextField(blank=True, null=True)
    source = models.TextField(blank=True, null=True)
    featuretype = models.TextField(blank=True, null=True)
    start = models.IntegerField(blank=True, null=True)
    end = models.IntegerField(blank=True, null=True)
    score = models.TextField(blank=True, null=True)
    strand = models.TextField(blank=True, null=True)
    frame = models.TextField(blank=True, null=True)
    attributes = models.TextField(blank=True, null=True)
    extra = models.TextField(blank=True, null=True)
    bin = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'features'


class Meta(models.Model):
    dialect = models.TextField(blank=True, null=True)
    version = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'meta'


class Relations(models.Model):
    parent = models.TextField(blank=True, null=True)
    child = models.TextField(blank=True, null=True)
    level = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'relations'


class SqliteStat1(models.Model):
    tbl = models.TextField(blank=True, null=True)  # This field type is a guess.
    idx = models.TextField(blank=True, null=True)  # This field type is a guess.
    stat = models.TextField(blank=True, null=True)  # This field type is a guess.

    class Meta:
        managed = False
        db_table = 'sqlite_stat1'
