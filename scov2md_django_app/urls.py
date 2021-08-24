from django.conf.urls import url,patterns,include
from django.conf import settings
from . import views

app_name = "covid19"

urlpatterns = [
    url(r'^search/$', views.index, name='index'),
    url(r'^project/$', views.project, name='project'),
    url(r'^$', views.home, name='home'),
    url(r'^home/$', views.home, name='home'),
    url(r'^prot/(?P<prot_name>\w+)/$', views.index, name='index'),
    url(r'^prot/(?P<prot_name>\w+)/(?P<genome_id>EPI_ISL_\d+)/$', views.index, name='index'),
    url(r'^allprot/(?P<genome_id>EPI_ISL_\d+)/$', views.index_allprot, name='index_allprot'),
    url(r'^(?P<dyn_id>[0-9]+)/$', views.dynanalysis, name='dynanalysis'),
    url(r'^(?P<dyn_id>[0-9]+)/(?P<sel_genome_id>EPI_ISL_\d+)/$', views.dynanalysis, name='dynanalysis'),
    url(r'^example/$', views.dynanalysis,{"variantimpact_def":True, "dyn_id":"28"}, name='dynanalysis'),
    url(r'^upload/$', views.upload, name='upload'),
    url(r'^upload/success/(?P<dyn_id>[0-9]+)$', views.upload_success, name='upload_success'),
    url(r'^ajax_notshow_warn/$', views.ajax_notshow_warn, name='ajax_notshow_warn'),
    url(r'^ajax_rmsd/$', views.ajax_rmsd, name='ajax_rmsd'),
    url(r'^ajax_rmsf/$', views.ajax_rmsf, name='ajax_rmsf'),
    url(r'^ajax_variant_impact/$', views.ajax_variant_impact, name='ajax_variant_impact'),
    url(r'^ajax_muts_in_isolate/$', views.ajax_muts_in_isolate, name='ajax_muts_in_isolate'),
    url(r'^ajax_autocomp_isolates/$', views.ajax_autocomp_isolates, name='ajax_autocomp_isolates'),
    url(r'^dwl/(?P<dyn_id>[0-9]+)/(?P<obj_type>[a-z]+)/(?P<obj_id>[0-9]+)/$', views.download_data, name="download_data"),
    url(r'^dwl/fasta/(?P<genome_id>EPI_ISL_\d+)/(?P<prot_name>\w+)/$', views.download_fasta, name="download_fasta"),
    url(r'^quickloadall/$', views.quickloadall, name="quickloadall"),
    url(r'^dwl/variantimpact/(?P<dyn_id>[0-9]+)/(?P<traj_id>[0-9]+)/(?P<position>[0-9]+)/(?P<analysis>\w+)/$', views.download_varimpact, name='download_varimpact'),
    url(r'^dwl/variantscores_traj/(?P<dyn_id>[0-9]+)/(?P<traj_id>[0-9]+)/(?P<protein>\w+)/(?P<position>[0-9]+)/(?P<variant>\w[0-9]+\w)/(?P<parameters_me>[\w_,]+)/(?P<parameters_td>[\w_,]+)/$', views.download_varscores_traj, name='download_varscores_traj'),
    url(r'^dwl/variantscores_all/(?P<dyn_id>[0-9]+)/$', views.download_variantscores_all, name='download_variantscores_all'),
#    url(r'^(?P<dyn_id>[0-9]+)/$', views.index, name='index'),
]

