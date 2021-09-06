from django.urls import path

from . import views

app_name = 'home'

urlpatterns = [
    path('about/', views.about),
    path('download/', views.download),
    path('', views.index, name='home'),
]