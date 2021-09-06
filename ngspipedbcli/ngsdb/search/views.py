from django.shortcuts import render

# Create your views here.
def index(request):
    return render(request, 'search/index.html')

# Create your views here.
def query(request):
    k = request.GET['keyword']
    k_dict = {'keywords': k}
    print(k_dict)
    return render(request, 'search/search_result_list.html', context=k_dict)