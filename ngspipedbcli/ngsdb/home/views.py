from django.shortcuts import render

# Create your views here.

def index(request):
    return render(request, 'home/index.html')

def about(request):
    return render(request, 'home/about.html')

def report(request):
    return render(request, 'home/report.html')

def download(request):
    return render(request, 'home/download.html')