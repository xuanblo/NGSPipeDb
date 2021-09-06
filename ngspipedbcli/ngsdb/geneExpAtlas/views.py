from django.shortcuts import render

import math
from django.http import JsonResponse, HttpResponse
from .models import Exp
from django.core import serializers
import json

# Create your views here.

def index(request):
    return render(request, 'geneExpAtlas/index.html')

def exp_json(request):
    # https://github.com/rg3915/django-datatables-experiment/issues/1
    # ?start=0&length=10
    exps = Exp.objects.all()
    total = exps.count()
    columns = [i.name for i in Exp._meta.get_fields()]

    _start = request.GET.get('start')
    _length = request.GET.get('length')
    if not _start:
        _start = '0'
    if not _length:
        _length = '10'
    if _start and _length:
        start = int(_start)
        length = int(_length)
        page = math.ceil(start / length) + 1 # math.ceil函数返回大于或等于一个给定数字的最小整数。
        per_page = length

        exps = exps[start: start+length]

    data = serializers.serialize("json", exps)
    data = json.loads(data)
    response = {
        'data': data,
        'page': page,  # [opcional]
        'per_page': per_page,  # [opcional]
        'recordsTotal': total,
        'recordsFiltered': total,
        'columns': columns,
        'columns_fields': [{'data':'pk'}]+[{'data':'fields.'+i} for i in columns[1:]],
    }
    
    return JsonResponse(response, content_type='application/json')

def exp_heatmap_json(request):
    import pandas as pd
    from clustergrammer import Network
    columns = [i.name for i in Exp._meta.get_fields()]
    #exps = Exp.objects.all().using("expDb").values_list("gene_id", "control_0", "control_1", "control_2", "treated_0", "treated_1", "treated_2")
    exps = Exp.objects.all().using("expDb").values()
    df = pd.DataFrame(list(exps), columns=columns)
    df.index = df.gene_id
    df = df.loc[:,df.columns[1:]]
    net = Network()
    net.load_df(df)

    # Z-score normalize the rows
    net.normalize(axis='row', norm_type='zscore', keep_orig=True)

    # filter for the top 100 columns based on their absolute value sum
    net.filter_N_top('col', 100, 'sum')

    # cluster using default parameters
    net.cluster()

    # save visualization JSON to file for use by front end
    data = net.export_net_json('viz')
    data = json.loads(data)
    #print(data)
    response = {
        'data': data,
    }
    return JsonResponse(response, content_type='application/json')

def heatmap(request):

    return render(request, 'geneExpAtlas/heatmap.html')