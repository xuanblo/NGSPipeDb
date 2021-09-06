from django.urls import path

from . import views

app_name = 'geneDetail'

urlpatterns = [
    path('', views.index),
]