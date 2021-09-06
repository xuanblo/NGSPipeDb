from django.urls import path

from . import views

app_name = 'geneExpAtlas'

urlpatterns = [
    path('', views.index, name='index'),
    path('heatmap', views.heatmap, name='heatmap'),
    path('exp_json', views.exp_json, name='exp_json'),
    path('exp_heatmap_json', views.exp_heatmap_json, name='exp_heatmap_json'),
]