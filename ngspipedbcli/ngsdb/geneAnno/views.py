from django.shortcuts import render

import math
from django.http import JsonResponse, HttpResponse
from .models import *
from django.core import serializers
import json

# Create your views here.

def index(request):
    return render(request, 'geneAnno/index.html')

def gff2json(request):
    # https://github.com/rg3915/django-datatables-experiment/issues/1
    # ?start=0&length=10
    features = Features.objects.all()
    total = features.count()
    
    _start = request.GET.get('start')
    _length = request.GET.get('length')
    if _start and _length:
        start = int(_start)
        length = int(_length)
        page = math.ceil(start / length) + 1 # math.ceil函数返回大于或等于一个给定数字的最小整数。
        per_page = length

        features = features[start: start+length]

    data = serializers.serialize("json", features)
    data = json.loads(data)
    response = {
        'data': data,
        'page': page,  # [opcional]
        'per_page': per_page,  # [opcional]
        'recordsTotal': total,
        'recordsFiltered': total,
    }
    
    return JsonResponse(response, content_type='application/json')

def anno_public_database(request):
    return render(request, 'geneAnno/public_database.html')

def anno_blast_genefamily(request):
    return render(request, 'geneAnno/blast2gene_family.html')