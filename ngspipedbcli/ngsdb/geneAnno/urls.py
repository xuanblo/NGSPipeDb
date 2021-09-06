from django.urls import path

from . import views

app_name = 'geneAnno'

urlpatterns = [
    path('', views.index, name='index'),
    path('gff2json', views.gff2json, name='gff2json'),
    path('anno_public_database', views.anno_public_database, name='anno_public_database'),
    path('anno_blast_genefamily', views.anno_blast_genefamily, name='anno_blast_genefamily'),
]