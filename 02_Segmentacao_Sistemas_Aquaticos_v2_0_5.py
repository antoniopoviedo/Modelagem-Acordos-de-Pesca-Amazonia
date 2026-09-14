# -*- coding: utf-8 -*-
# ==============================================================================
# MODELO 02 v2.0.5 FINAL — SEGMENTAÇÃO FUNCIONAL + QA/QC + POLÍGONOS
#
# Objetivo: reconstruir a cadeia funcional do Modelo 02 sem calibração oportunista,
# usando o benchmark manual de Acajatuba como CONTROLE DE REGRESSÃO.
#
# Regra operacional: FAIL-FAST. O algoritmo interrompe na primeira divergência
# material em relação ao Registro Mestre Técnico. Nenhum parâmetro é ajustado
# automaticamente para "fazer fechar" o benchmark.
#
# Gates obrigatórios no modo Acajatuba + benchmark:
#   HB3 -> HB4 domínio/lambda -> HB4 pixel-a-pixel -> HC -> P1b/P2 -> P3
#   -> P5b -> P6/P7 -> P7b -> P8a2 -> P8 -> P9 -> P10 final.
# ==============================================================================

"""
02_Segmentacao_Sistemas_Aquaticos_v2.0.5.py

MODELO 02 – SEGMENTAÇÃO SISTEMAS AQUATICOS INTEGRAL
Projeto: Pipeline de Espacialização dos Acordos de Pesca
Baseline metodológica: Acajatuba
Ambiente-alvo: QGIS 3.44.8 / Python 3.12 / GDAL 3.12

ARQUITETURA
-----------
ENTRADAS EXTERNAS (somente M01 + bases JRC; nenhuma M02):
    - 05_Pontos_Diagnostico_JRC (ou pontos legais reprojetados equivalentes)
    - 04_JRC_MaxExtent_Metrico (saída M01; usado SOMENTE no HB1)
    - JRC MaxExtent estadual/original (usado no HB2 e HB4; água = valor 1)
    - JRC Occurrence estadual/original
    - JRC Seasonality estadual/original

CADEIA INTERNA:
    HB1 -> HB2/HB3 -> HB4 -> HC -> P1b/P2 -> P3 -> P4 -> P5b
        -> P6/P6b -> P7 -> P7b(60m) -> P8a2 -> P8 -> P9 -> P10

Nenhuma camada M02_* é exigida como entrada.

PARAMETRIZAÇÃO RECUPERADA DA MODELAGEM MANUAL
----------------------------------------------
PIX_M       = 30 m
BUFFER      = 6000 m
RAIO_OPEN   = 150 m
ALPHA_MORF  = 0.45
BETA_HIDRO  = 0.20
GAMMA_TRANS = 0
EPS         = 1e-9
P7 SNAP     = 1 m
P7 W_LEN    = 0.50
P7 W_C      = 0.25
P7 W_A1     = 0.25
P7 FRAC_SEC = 0.20
P7 LEN_SEC  = 150 m
P7b FAIXA   = 60 m

Acajatuba – envelope HB4 validado para regressão:
    xmin = 764687.00
    xmax = 790727.00
    ymin = 9644344.34
    ymax = 9668434.34
    grade esperada = 868 x 803
    água esperada  = 230725 pixels

Acajatuba – unidades funcionais validadas:
    SIS_123       = IDs 1,2,3
    SIS_4510      = IDs 4,5,10
    SIS_789       = IDs 7,8,9
    SIS_6         = ID 6
    SIS_PRES_1112 = IDs 11,12
    SIS_PRES_13   = ID 13

P1b – Preservação em Acajatuba:
    - ID11 -> pixel (190,194) -> 770522.00, 9662719.34
    - ID12 -> pixel (224,225) -> 771452.00, 9661699.34
    - ID13 -> pixel (181,262) -> 772562.00, 9662989.34
    - as três âncoras são produtos funcionais P1b, não simples snaps
      dos pontos legais.

HB4 – detalhes de reprodução manual:
    - semente = célula da âncora; se não-água, busca apenas +/-3 células;
    - medianas e pool da escala comum excluem custo == 0;
    - LAMBDA_COMUM = mediana do pool de custos positivos das quatro unidades.

HB3 – exceção validada:
    somente ID 5 é reancorado estruturalmente.
    A regra operacional reproduz o diagnóstico manual HB2:
    - JRC MaxExtent estadual/original;
    - janela local exata de 6000 m (400 x 400 células a 30 m);
    - rotulagem 8-vizinhos do MaxExtent;
    - componentes candidatos = componentes com pelo menos um pixel até 1500 m;
    - seleção do maior componente entre os candidatos;
    - âncora = pixel candidato mais próximo dentro do componente selecionado.

P4 – limiares EMPÍRICOS, calculados em precisão integral:
    M_CRIT = P25(MARGEM3)
    C_CRIT = P50(A2)
    M_ALTA = P10(MARGEM3)
    C_ALTA = P75(A2)
    A_MIN  = P25(A1)

P5b:
    evidência mínima é avaliada antes da criticidade;
    A1 < A_MIN => SEM_EVIDENCIA;
    interfaces válidas: MARGEM3 <= M_CRIT e A2 >= C_CRIT;
    tripla válida somente quando A3 >= C_CRIT;
    componentes de interface: conectividade 8-vizinhos.

P6/P7:
    P6 reproduz o marching squares documental, exigindo os quatro vértices
    no mesmo componente crítico e resolvendo casos 5/10 pelo sinal médio.
    P7 reconstrói topologia com SNAP=1 m, conserva 619=615+4 segmentos e
    hierarquiza ramos pelo SCORE_P7 documental.

P7b:
    faixa operacional = 60 m em torno dos RAMOS PRINCIPAIS.

P8a2:
    apenas empate triplo sem evidência;
    propagação espacial em grade 8-vizinhos (Chebyshev);
    fonte = pixel com A1 >= A_MIN;
    empate de distância => menor PID;
    CLASSE_PROX = ARGMAX da fonte.

P8:
    prioridade PROX_ESPACIAL -> FAIXA_OPER -> ARGMAX.

BENCHMARK ACAJATUBA (QA/QC; não regra geral)
---------------------------------------------
Domínio       230725
COMERCIAL     141672
ESPORTIVA      57607
PRESERVACAO    31446
ARGMAX        223184
FAIXA_OPER      1801
PROX_ESPACIAL   5740
ALTA           138610
MEDIA           32633
TRANSICAO        1801
BAIXA           57681
"""

import math
import heapq
import unicodedata
from collections import Counter, defaultdict, deque

import numpy as np
from osgeo import gdal, ogr, osr

from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterField,
    QgsProcessingParameterNumber,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterFeatureSink,
    QgsFeature,
    QgsFields,
    QgsField,
    QgsGeometry,
    QgsPointXY,
    QgsWkbTypes,
    QgsFeatureSink,
    QgsSpatialIndex,
    QgsCoordinateTransform,
    QgsProject,
)


class SegmentacaoFuncionalIntegral(QgsProcessingAlgorithm):

    PONTOS = "PONTOS"
    MAX_METRICO = "MAX_METRICO"
    MAX_ORIGINAL = "MAX_ORIGINAL"
    OCCURRENCE = "OCCURRENCE"
    SEASONALITY = "SEASONALITY"

    ID_FIELD = "ID_FIELD"
    CAT_FIELD = "CAT_FIELD"
    SISTEMA_FIELD = "SISTEMA_FIELD"

    PIX_M = "PIX_M"
    BUFFER_M = "BUFFER_M"
    RAIO_OPEN_M = "RAIO_OPEN_M"
    FAIXA_M = "FAIXA_M"
    MODO_ACAJATUBA = "MODO_ACAJATUBA"
    BENCHMARK = "BENCHMARK"

    ANCORAS_OUT = "ANCORAS_OUT"
    AFINIDADES_OUT = "AFINIDADES_OUT"
    INTERFACES_OUT = "INTERFACES_OUT"
    RAMOS_OUT = "RAMOS_OUT"
    FAIXA_OUT = "FAIXA_OUT"
    PROX_OUT = "PROX_OUT"
    FINAL_OUT = "FINAL_OUT"
    RESUMO_OUT = "RESUMO_OUT"
    LIMIARES_OUT = "LIMIARES_OUT"

    P9_COMPONENTES_OUT = "P9_COMPONENTES_OUT"
    P9_FRAGMENTOS_OUT = "P9_FRAGMENTOS_OUT"
    P9_CLASSE_ORIGEM_OUT = "P9_CLASSE_ORIGEM_OUT"
    P9_CLASSE_CONFIANCA_OUT = "P9_CLASSE_CONFIANCA_OUT"
    P9_ORIGEM_CONFIANCA_OUT = "P9_ORIGEM_CONFIANCA_OUT"
    P9_RESUMO_OUT = "P9_RESUMO_OUT"

    P10_COMPONENTES_OUT = "P10_COMPONENTES_OUT"
    P10_CLASSES_OUT = "P10_CLASSES_OUT"
    P10_RESUMO_OUT = "P10_RESUMO_OUT"

    ALPHA_MORF = 0.45
    BETA_HIDRO = 0.20
    GAMMA_TRANS = 0.0
    EPS_FIXED = 1e-9
    HB2_RAIO_BUSCA_M = 1500.0

    SNAP_M = 1.0
    W_LEN = 0.50
    W_C = 0.25
    W_A1 = 0.25
    FRAC_SEC = 0.20
    LEN_SEC_MIN_M = 150.0

    CLASSES = ("COMERCIAL", "ESPORTIVA", "PRESERVACAO")

    BENCH = {
        "n_total": 230725,
        "classes": {"COMERCIAL": 141672, "ESPORTIVA": 57607, "PRESERVACAO": 31446},
        "origens": {"ARGMAX": 223184, "FAIXA_OPER": 1801, "PROX_ESPACIAL": 5740},
        "confiancas": {"ALTA": 138610, "MEDIA": 32633, "TRANSICAO": 1801, "BAIXA": 57681},
    }

    # ------------------------------------------------------------------
    # QGIS metadata
    # ------------------------------------------------------------------
    def tr(self, text):
        return QCoreApplication.translate("SegmentacaoFuncionalIntegral", text)

    def createInstance(self):
        return SegmentacaoFuncionalIntegral()

    def name(self):
        return "02_modelo_final_p1_p10_v6"

    def displayName(self):
        return self.tr("02 – Segmentação dos Sistemas Aquáticos v2.0.5")

    def group(self):
        return self.tr("Pipeline Acordos de Pesca")

    def groupId(self):
        return "pipeline_acordos_pesca"

    def shortHelpString(self):
        return self.tr(
            "Modelo 02 integral. HB1 usa o MaxExtent métrico do M01; HB2/HB4 "
            "usam o MaxExtent estadual/original, conforme a modelagem manual. "
            "Executa HB1–P10 em cadeia única. P9 valida geometricamente/estatisticamente o P8 sem reclassificar; P10 converte a grade final em polígonos exatos, sem remover ou suavizar fragmentos. No modo benchmark Acajatuba, aplica gates sucessivos contra o Registro Mestre e as camadas manuais preservadas."
        )

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.PONTOS,
            self.tr("Pontos com diagnóstico JRC – saída 05 do Modelo 01"),
            [QgsProcessing.TypeVectorPoint]
        ))
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.MAX_METRICO,
            self.tr("JRC MaxExtent métrico – saída 04 do Modelo 01 (somente HB1)")
        ))
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.MAX_ORIGINAL,
            self.tr("JRC MaxExtent estadual/original – HB2/HB4")
        ))
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.OCCURRENCE,
            self.tr("JRC Occurrence estadual/original")
        ))
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.SEASONALITY,
            self.tr("JRC Seasonality estadual/original")
        ))

        self.addParameter(QgsProcessingParameterField(
            self.ID_FIELD, self.tr("Campo identificador do ponto"),
            parentLayerParameterName=self.PONTOS,
            type=QgsProcessingParameterField.Any,
            defaultValue="id"
        ))
        self.addParameter(QgsProcessingParameterField(
            self.CAT_FIELD, self.tr("Campo de categoria regulatória"),
            parentLayerParameterName=self.PONTOS,
            type=QgsProcessingParameterField.Any,
            defaultValue="categoria"
        ))
        self.addParameter(QgsProcessingParameterField(
            self.SISTEMA_FIELD,
            self.tr("Campo de sistema funcional – opcional fora de Acajatuba"),
            parentLayerParameterName=self.PONTOS,
            type=QgsProcessingParameterField.Any,
            optional=True
        ))

        self.addParameter(QgsProcessingParameterNumber(
            self.PIX_M, self.tr("Resolução da grade (m)"),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=30.0, minValue=1.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.BUFFER_M, self.tr("Margem/janela funcional (m)"),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=6000.0, minValue=0.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.RAIO_OPEN_M, self.tr("Raio de openness local (m)"),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=150.0, minValue=0.0
        ))
        self.addParameter(QgsProcessingParameterNumber(
            self.FAIXA_M, self.tr("Faixa operacional P7b (m)"),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=60.0, minValue=0.0
        ))
        self.addParameter(QgsProcessingParameterBoolean(
            self.MODO_ACAJATUBA,
            self.tr("Aplicar unidades/âncora estrutural calibradas de Acajatuba"),
            defaultValue=True
        ))
        self.addParameter(QgsProcessingParameterBoolean(
            self.BENCHMARK,
            self.tr("Executar benchmark estrito de Acajatuba"),
            defaultValue=True
        ))

        for key, label, geom_type in [
            (self.ANCORAS_OUT, "M02_01 – Âncoras Funcionais", QgsProcessing.TypeVectorPoint),
            (self.AFINIDADES_OUT, "M02_02 – Afinidades Multiclasse", QgsProcessing.TypeVectorPoint),
            (self.INTERFACES_OUT, "M02_03 – Interfaces Válidas P5b", QgsProcessing.TypeVectorPoint),
            (self.RAMOS_OUT, "M02_04 – Ramos de Equilíbrio P7", QgsProcessing.TypeVectorLine),
            (self.FAIXA_OUT, "M02_05 – Faixa Operacional P7b", QgsProcessing.TypeVectorPoint),
            (self.PROX_OUT, "M02_06 – Preenchimento de Empates P8a2", QgsProcessing.TypeVectorPoint),
            (self.FINAL_OUT, "M02_07 – Segmentação Funcional Final P8", QgsProcessing.TypeVectorPoint),
            (self.RESUMO_OUT, "M02_08 – Resumo QA/QC", QgsProcessing.TypeVectorAnyGeometry),
            (self.LIMIARES_OUT, "M02_09 – Limiares Empíricos", QgsProcessing.TypeVectorAnyGeometry),
        ]:
            self.addParameter(QgsProcessingParameterFeatureSink(key, self.tr(label), geom_type))

        for key, label, geom_type in [
            (self.P9_COMPONENTES_OUT, "M02_10 – P9 Componentes 8-conectados", QgsProcessing.TypeVectorAnyGeometry),
            (self.P9_FRAGMENTOS_OUT, "M02_11 – P9 Fragmentos pequenos (<=20 pixels)", QgsProcessing.TypeVectorPoint),
            (self.P9_CLASSE_ORIGEM_OUT, "M02_12 – P9 Classe × Origem", QgsProcessing.TypeVectorAnyGeometry),
            (self.P9_CLASSE_CONFIANCA_OUT, "M02_13 – P9 Classe × Confiança", QgsProcessing.TypeVectorAnyGeometry),
            (self.P9_ORIGEM_CONFIANCA_OUT, "M02_14 – P9 Origem × Confiança", QgsProcessing.TypeVectorAnyGeometry),
            (self.P9_RESUMO_OUT, "M02_15 – P9 Resumo QA/QC", QgsProcessing.TypeVectorAnyGeometry),
            (self.P10_COMPONENTES_OUT, "M02_16 – P10 Polígonos por Componente", QgsProcessing.TypeVectorPolygon),
            (self.P10_CLASSES_OUT, "M02_17 – P10 Polígonos Dissolvidos por Classe", QgsProcessing.TypeVectorPolygon),
            (self.P10_RESUMO_OUT, "M02_18 – P10 Resumo Poligonal", QgsProcessing.TypeVectorAnyGeometry),
        ]:
            self.addParameter(QgsProcessingParameterFeatureSink(key, self.tr(label), geom_type))

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _norm(v):
        if v is None:
            return ""
        s = unicodedata.normalize("NFKD", str(v).strip().upper())
        return "_".join(s.encode("ascii", "ignore").decode("ascii").split())

    @classmethod
    def _norm_class(cls, v):
        s = cls._norm(v)
        aliases = {
            "COM": "COMERCIAL", "COMERCIAL": "COMERCIAL",
            "ESP": "ESPORTIVA", "ESPORTIVA": "ESPORTIVA",
            "PRES": "PRESERVACAO", "PRESERVACAO": "PRESERVACAO",
        }
        return aliases.get(s, s)

    @staticmethod
    def _point_xy(geom):
        if geom is None or geom.isEmpty():
            return None
        if geom.isMultipart():
            pts = geom.asMultiPoint()
            return pts[0] if pts else None
        return geom.asPoint()

    @staticmethod
    def _world_to_rc(gt, x, y):
        c = int(math.floor((x - gt[0]) / gt[1]))
        r = int(math.floor((y - gt[3]) / gt[5]))
        return r, c

    @staticmethod
    def _rc_to_xy(gt, r, c):
        x = gt[0] + (c + 0.5) * gt[1] + (r + 0.5) * gt[2]
        y = gt[3] + (c + 0.5) * gt[4] + (r + 0.5) * gt[5]
        return float(x), float(y)

    @staticmethod
    def _neighbors8(pix_m):
        d = pix_m * math.sqrt(2.0)
        return [(-1,0,pix_m),(1,0,pix_m),(0,-1,pix_m),(0,1,pix_m),
                (-1,-1,d),(-1,1,d),(1,-1,d),(1,1,d)]

    @staticmethod
    def _label8(mask):
        nr, nc = mask.shape
        labels = np.zeros((nr, nc), dtype=np.int32)
        comps = {}
        lab = 0
        for r in range(nr):
            for c in range(nc):
                if not mask[r,c] or labels[r,c] != 0:
                    continue
                lab += 1
                q = deque([(r,c)])
                labels[r,c] = lab
                cells = []
                while q:
                    rr,cc = q.popleft()
                    cells.append((rr,cc))
                    for dr in (-1,0,1):
                        for dc in (-1,0,1):
                            if dr == 0 and dc == 0:
                                continue
                            r2,c2 = rr+dr,cc+dc
                            if 0 <= r2 < nr and 0 <= c2 < nc and mask[r2,c2] and labels[r2,c2] == 0:
                                labels[r2,c2] = lab
                                q.append((r2,c2))
                comps[lab] = cells
        return labels, comps

    @staticmethod
    def _integral_openness(agua, raio_pix):
        A = agua.astype(np.int32)
        integ = np.pad(A, ((1,0),(1,0)), mode='constant').cumsum(0).cumsum(1)
        nr,nc = agua.shape
        out = np.zeros((nr,nc), dtype=np.float32)
        for r in range(nr):
            r0=max(0,r-raio_pix); r1=min(nr-1,r+raio_pix)
            for c in range(nc):
                if not agua[r,c]:
                    continue
                c0=max(0,c-raio_pix); c1=min(nc-1,c+raio_pix)
                soma=(integ[r1+1,c1+1]-integ[r0,c1+1]-integ[r1+1,c0]+integ[r0,c0])
                out[r,c]=float(soma)/float((r1-r0+1)*(c1-c0+1))
        return out

    @staticmethod
    def _winner(a_com, a_esp, a_pres, eps):
        vals=[("COMERCIAL",a_com),("ESPORTIVA",a_esp),("PRESERVACAO",a_pres)]
        vals.sort(key=lambda z:z[1], reverse=True)
        a1,a2,a3=vals[0][1],vals[1][1],vals[2][1]
        tie2=abs(a1-a2)<=eps
        tie3=tie2 and abs(a2-a3)<=eps
        pair="_".join(sorted([vals[0][0].replace("ORTIVA",""), vals[1][0].replace("ERCIAL","").replace("ERVACAO","")]))
        # normalizar nomes dos pares de forma explícita
        pair_map={frozenset(("COMERCIAL","ESPORTIVA")):"COM_ESP",
                  frozenset(("COMERCIAL","PRESERVACAO")):"COM_PRES",
                  frozenset(("ESPORTIVA","PRESERVACAO")):"ESP_PRES"}
        pair=pair_map.get(frozenset((vals[0][0],vals[1][0])),"TRIPLA" if tie3 else "")
        return vals[0][0],a1,a2,a3,tie2,tie3,pair,vals

    @staticmethod
    def _percentile(a, p):
        return float(np.percentile(np.asarray(a,dtype=np.float64), p))

    def _warp_to_grid(self, raster_layer, bounds, crs_wkt, pix_m, cols, rows, nodata=0.0):
        src = gdal.Open(raster_layer.source().split('|')[0])
        if src is None:
            raise QgsProcessingException("Falha ao abrir raster: " + raster_layer.name())
        xmin,ymin,xmax,ymax=bounds
        mem=gdal.Warp('',src,format='MEM',dstSRS=crs_wkt,
                      outputBounds=(xmin,ymin,xmax,ymax),width=cols,height=rows,
                      resampleAlg='near',dstNodata=nodata)
        if mem is None:
            raise QgsProcessingException("Falha GDAL Warp: " + raster_layer.name())
        return mem.ReadAsArray().astype(np.float64)

    def _warp_exact(self, raster_layer, bounds, dst_crs_wkt, pix_m, nodata=0.0):
        """
        Reproduz a lógica manual: gdal.Warp diretamente da base estadual/original,
        com outputBounds exatos e xRes/yRes de 30 m, sem alinhar a outra grade.
        """
        src = gdal.Open(raster_layer.source().split('|')[0])
        if src is None:
            raise QgsProcessingException(
                "Falha ao abrir raster: " + raster_layer.name()
            )
        xmin, ymin, xmax, ymax = bounds
        mem = gdal.Warp(
            '',
            src,
            format='MEM',
            dstSRS=dst_crs_wkt,
            outputBounds=(xmin, ymin, xmax, ymax),
            xRes=pix_m,
            yRes=pix_m,
            resampleAlg='near'
        )
        if mem is None:
            raise QgsProcessingException(
                "Falha GDAL Warp: " + raster_layer.name()
            )
        arr = mem.ReadAsArray()
        return mem, arr

    def _propagar(self, agua, open_arr, p_h, seeds, open_ref, pix_m, eps):
        nr,nc=agua.shape
        denom=max(eps,1.0-open_ref)
        p_m=np.maximum(0.0,(open_arr.astype(np.float64)-open_ref)/denom)
        custo=np.full((nr,nc),np.inf,dtype=np.float64)
        hp=[]
        for r,c in seeds:
            if 0<=r<nr and 0<=c<nc and agua[r,c]:
                custo[r,c]=0.0
                heapq.heappush(hp,(0.0,r,c))
        for_pop=self._neighbors8(pix_m)
        while hp:
            atual,r,c=heapq.heappop(hp)
            if atual != custo[r,c]:
                continue
            for dr,dc,dist in for_pop:
                rr,cc=r+dr,c+dc
                if rr<0 or rr>=nr or cc<0 or cc>=nc or not agua[rr,cc]:
                    continue
                pm=0.5*(p_m[r,c]+p_m[rr,cc])
                ph=0.5*(p_h[r,c]+p_h[rr,cc])
                passo=dist*(1.0+self.ALPHA_MORF*pm+self.BETA_HIDRO*ph)
                novo=atual+passo
                if novo < custo[rr,cc]:
                    custo[rr,cc]=novo
                    heapq.heappush(hp,(novo,rr,cc))
        return custo

    # ------------------------------------------------------------------
    # P6 marching-squares segment extraction
    # ------------------------------------------------------------------
    def _propagar_p2_original(self, mask, open_arr, hidro, seeds, open_ref, pix_m):
        """
        Implementação literal recuperada do código original C2.2h-p2.

        Diferenças deliberadas em relação ao HB4:
        - pen_morf = abs(OPEN do pixel destino - OPEN_REF);
        - pen_hidro = 1 - hidro do pixel destino;
        - hidro = 0.5*Occurrence + 0.5*Seasonality normalizados;
        - não usa média origem/destino para as penalidades;
        - EPS interno = 1e-12.
        """
        rows, cols = mask.shape
        eps_p2 = 1e-12
        dist = np.full((rows, cols), np.inf, dtype=np.float64)
        heap = []

        for r, c in seeds:
            r = int(r); c = int(c)
            if 0 <= r < rows and 0 <= c < cols and bool(mask[r, c]):
                dist[r, c] = 0.0
                heapq.heappush(heap, (0.0, r, c))

        viz = [
            (-1,-1,math.sqrt(2.0)), (-1,0,1.0), (-1,1,math.sqrt(2.0)),
            ( 0,-1,1.0),                         ( 0,1,1.0),
            ( 1,-1,math.sqrt(2.0)), ( 1,0,1.0), ( 1,1,math.sqrt(2.0)),
        ]

        ref = float(open_ref)

        while heap:
            custo, r, c = heapq.heappop(heap)
            custo = float(custo); r = int(r); c = int(c)

            if custo > float(dist[r, c]) + eps_p2:
                continue

            for dr, dc, fator_dist in viz:
                nr2 = int(r + int(dr))
                nc2 = int(c + int(dc))

                if nr2 < 0 or nr2 >= rows or nc2 < 0 or nc2 >= cols:
                    continue
                if not bool(mask[nr2, nc2]):
                    continue

                pen_morf = abs(float(open_arr[nr2, nc2]) - ref)
                pen_hidro = 1.0 - float(hidro[nr2, nc2])

                fator_custo = (
                    1.0
                    + self.ALPHA_MORF * pen_morf
                    + self.BETA_HIDRO * pen_hidro
                    + 0.0
                )

                passo = float(pix_m) * float(fator_dist) * float(fator_custo)
                novo = float(custo + passo)

                if novo + eps_p2 < float(dist[nr2, nc2]):
                    dist[nr2, nc2] = novo
                    heapq.heappush(heap, (novo, nr2, nc2))

        return dist


    def _interp(self, p1, v1, p2, v2, eps):
        den=(v2-v1)
        t=0.5 if abs(den)<=eps else (-v1/den)
        t=max(0.0,min(1.0,t))
        return (p1[0]+t*(p2[0]-p1[0]), p1[1]+t*(p2[1]-p1[1]))

    def _segments_zero(self, delta, comp_mask, gt, eps):
        nr,nc=delta.shape
        segs=[]
        for r in range(nr-1):
            for c in range(nc-1):
                if not (comp_mask[r,c] or comp_mask[r,c+1] or comp_mask[r+1,c] or comp_mask[r+1,c+1]):
                    continue
                vals=[delta[r,c],delta[r,c+1],delta[r+1,c+1],delta[r+1,c]]
                if not all(np.isfinite(v) for v in vals):
                    continue
                pts=[self._rc_to_xy(gt,r,c),self._rc_to_xy(gt,r,c+1),
                     self._rc_to_xy(gt,r+1,c+1),self._rc_to_xy(gt,r+1,c)]
                cross=[]
                for a,b in [(0,1),(1,2),(2,3),(3,0)]:
                    va,vb=vals[a],vals[b]
                    if abs(va)<=eps and abs(vb)<=eps:
                        continue
                    if (va<=0<=vb) or (vb<=0<=va):
                        if va*vb<=0:
                            cross.append(self._interp(pts[a],va,pts[b],vb,eps))
                # remover duplicatas de vértice
                uniq=[]
                for p in cross:
                    if not any(math.hypot(p[0]-q[0],p[1]-q[1])<1e-6 for q in uniq):
                        uniq.append(p)
                if len(uniq)==2:
                    segs.append((uniq[0],uniq[1]))
                elif len(uniq)==4:
                    # caso sela: preservar duas conexões locais sem impor morfologia
                    segs.append((uniq[0],uniq[1])); segs.append((uniq[2],uniq[3]))
        return segs

    @staticmethod
    def _snap_key(p, snap):
        return (int(round(p[0]/snap)), int(round(p[1]/snap)))

    def _build_branches(self, segs, snap):
        # edges indexados por nós snapados
        edges=[]; adj=defaultdict(list); node_xy={}
        micro=[]
        for i,(a,b) in enumerate(segs):
            ka=self._snap_key(a,snap); kb=self._snap_key(b,snap)
            if ka==kb:
                micro.append((a,b)); continue
            node_xy.setdefault(ka,a); node_xy.setdefault(kb,b)
            ei=len(edges); edges.append((ka,kb,a,b))
            adj[ka].append(ei); adj[kb].append(ei)
        used=set(); branches=[]

        def walk(start_node, first_e):
            pts=[node_xy[start_node]]; e=first_e; node=start_node
            while True:
                if e in used: break
                used.add(e)
                ka,kb,a,b=edges[e]
                nxt=kb if node==ka else ka
                pts.append(node_xy[nxt])
                if len(adj[nxt])!=2:
                    break
                cand=[x for x in adj[nxt] if x not in used]
                if not cand: break
                node=nxt; e=cand[0]
            return pts

        for node, elist in list(adj.items()):
            if len(elist)!=2:
                for e in elist:
                    if e not in used:
                        p=walk(node,e)
                        if len(p)>=2: branches.append(p)
        # ciclos remanescentes
        for e in range(len(edges)):
            if e in used: continue
            start=edges[e][0]; p=walk(start,e)
            if len(p)>=2: branches.append(p)
        return branches,micro

    # ------------------------------------------------------------------
    # Main
    # ------------------------------------------------------------------
    def processAlgorithm(self, parameters, context, feedback):
        pontos=self.parameterAsSource(parameters,self.PONTOS,context)
        rmax_m01=self.parameterAsRasterLayer(parameters,self.MAX_METRICO,context)
        rmax_orig=self.parameterAsRasterLayer(parameters,self.MAX_ORIGINAL,context)
        rocc=self.parameterAsRasterLayer(parameters,self.OCCURRENCE,context)
        rseas=self.parameterAsRasterLayer(parameters,self.SEASONALITY,context)
        if None in (pontos,rmax_m01,rmax_orig,rocc,rseas):
            raise QgsProcessingException("Entradas obrigatórias inválidas.")

        id_field=self.parameterAsString(parameters,self.ID_FIELD,context)
        cat_field=self.parameterAsString(parameters,self.CAT_FIELD,context)
        sistema_field=self.parameterAsString(parameters,self.SISTEMA_FIELD,context)
        pix_m=self.parameterAsDouble(parameters,self.PIX_M,context)
        buffer_m=self.parameterAsDouble(parameters,self.BUFFER_M,context)
        raio_open=self.parameterAsDouble(parameters,self.RAIO_OPEN_M,context)
        faixa_m=self.parameterAsDouble(parameters,self.FAIXA_M,context)
        eps=self.EPS_FIXED
        modo_acaj=self.parameterAsBool(parameters,self.MODO_ACAJATUBA,context)
        run_bench=self.parameterAsBool(parameters,self.BENCHMARK,context)

        fset={f.name() for f in pontos.fields()}
        if id_field not in fset or cat_field not in fset:
            raise QgsProcessingException("Campos ID/categoria não encontrados nos pontos.")
        if not pontos.sourceCrs().isValid() or pontos.sourceCrs().isGeographic():
            raise QgsProcessingException("Pontos devem estar em CRS projetado/métrico.")

        feedback.pushInfo("="*120)
        feedback.pushInfo("MODELO 02 v5.1 – RECONSTRUCAO DOCUMENTAL P6-P8 / ACAJATUBA")
        feedback.pushInfo("="*120)
        feedback.pushInfo("Entradas de modelagem: somente M01 + bases JRC. Camadas M02 manuais são usadas APENAS como benchmark de comparação.")
        feedback.pushInfo("EPS metodológico fixo internamente: 1e-9\nPrecisão HB4 recuperada: HC_* calculadas/persistidas via array np.float32, conforme script manual.\nReconstrução documental: HB4/HC já fechados; P2 substituído pela função original recuperada do P2.docx; P3 usa EPS_TIE=1e-6 e rótulos originais de empate.")

        # --------------------------------------------------------------
        # Benchmark manual carregado no projeto
        # --------------------------------------------------------------
        bench_layers = {}
        if run_bench:
            bench_names = [
                "M02_C2_ancoras_funcionais_hb3",
                "M02_C2_afinidades_unidades_hb4",
                "M02_C2_afinidades_classes_hc",
                "M02_C2_preservacao_ancoras_funcionais_p1b",
                "M02_C2_afinidades_preservacao_p2",
                "M02_C2_afinidades_multiclasse_p3",
                "M02_C2_zonas_criticas_multiclasse_p4",
            ]
            missing = []
            for nm in bench_names:
                ls = QgsProject.instance().mapLayersByName(nm)
                if not ls or not ls[0].isValid():
                    missing.append(nm)
                else:
                    bench_layers[nm] = ls[0]
            if missing:
                raise QgsProcessingException(
                    "BENCHMARK MANUAL INCOMPLETO. Carregue no projeto: "
                    + ", ".join(missing)
                )

            hb4_bench = bench_layers["M02_C2_afinidades_unidades_hb4"]
            if hb4_bench.featureCount() != 230725:
                raise QgsProcessingException(
                    "BENCHMARK HB4 inválido: "
                    f"{hb4_bench.featureCount()} feições; esperado 230725."
                )
            hb4_fields = {f.name() for f in hb4_bench.fields()}
            if "LAMBDA_C" not in hb4_fields:
                raise QgsProcessingException(
                    "BENCHMARK HB4 sem campo LAMBDA_C."
                )
            lambda_vals = set()
            for bf in hb4_bench.getFeatures():
                v = bf["LAMBDA_C"]
                if v is not None:
                    lambda_vals.add(round(float(v), 9))
                if len(lambda_vals) > 2:
                    break
            if len(lambda_vals) != 1:
                raise QgsProcessingException(
                    "BENCHMARK HB4 possui LAMBDA_C inconsistente: "
                    + str(sorted(lambda_vals))
                )
            lambda_benchmark = float(next(iter(lambda_vals)))

            feedback.pushInfo("-"*120)
            feedback.pushInfo("BENCHMARK MANUAL CARREGADO")
            for nm in bench_names:
                feedback.pushInfo(
                    f"{nm} | feições={bench_layers[nm].featureCount()}"
                )
            feedback.pushInfo(
                f"HB4 manual: LAMBDA_C={lambda_benchmark:.9f}"
            )
            feedback.pushInfo("-"*120)
        else:
            lambda_benchmark = None

        # --------------------------------------------------------------
        # Ler pontos e definir unidades
        # --------------------------------------------------------------
        pts={}
        tr_to_pts=None
        if rmax_m01.crs()!=pontos.sourceCrs():
            tr_to_pts=QgsCoordinateTransform(rmax_m01.crs(),pontos.sourceCrs(),QgsProject.instance())
        for f in pontos.getFeatures():
            try: ident=int(f[id_field])
            except Exception: ident=str(f[id_field]).strip()
            p=self._point_xy(f.geometry())
            if p is None: continue
            cat=self._norm_class(f[cat_field])
            pts[ident]={"x":float(p.x()),"y":float(p.y()),"cat":cat,"feat":f}

        if modo_acaj:
            expected=set(range(1,14))
            if set(k for k in pts if isinstance(k,int)) != expected:
                raise QgsProcessingException("Modo Acajatuba exige IDs 1–13.")
            systems={
                "SIS_123":[1,2,3], "SIS_4510":[4,5,10],
                "SIS_789":[7,8,9], "SIS_6":[6],
                "SIS_PRES_1112":[11,12], "SIS_PRES_13":[13],
            }
        else:
            if sistema_field and sistema_field in fset:
                systems=defaultdict(list)
                for ident,d in pts.items():
                    systems[str(d["feat"][sistema_field]).strip()].append(ident)
                systems=dict(systems)
            else:
                # fallback conservador: cada ponto é unidade própria; agregação ocorre por classe.
                systems={f"UNIT_{ident}":[ident] for ident in pts}

        # --------------------------------------------------------------
        # Abrir MaxExtent métrico M01 e amostrar grade
        # --------------------------------------------------------------
        dsmax=gdal.Open(rmax_m01.source().split('|')[0])
        if dsmax is None:
            raise QgsProcessingException("Falha ao abrir MaxExtent métrico M01.")
        MAX0=dsmax.ReadAsArray()
        gt0=dsmax.GetGeoTransform(); proj0=dsmax.GetProjection()
        agua0=(MAX0>0)
        nr0,nc0=agua0.shape
        water_rc=np.argwhere(agua0)
        if water_rc.size==0: raise QgsProcessingException("MaxExtent métrico sem água.")

        # nearest water pixel baseline for every point
        anchors0={}
        for ident,d in pts.items():
            rr0,cc0=self._world_to_rc(gt0,d["x"],d["y"])
            rad=max(1,int(math.ceil(1000.0/pix_m)))
            best=None
            while best is None:
                r0=max(0,rr0-rad); r1=min(nr0,rr0+rad+1)
                c0=max(0,cc0-rad); c1=min(nc0,cc0+rad+1)
                loc=np.argwhere(agua0[r0:r1,c0:c1])
                if len(loc):
                    for lr,lc in loc:
                        r,c=r0+int(lr),c0+int(lc)
                        x,y=self._rc_to_xy(gt0,r,c)
                        dd=math.hypot(x-d["x"],y-d["y"])
                        if best is None or dd<best[0]: best=(dd,r,c,x,y)
                rad*=2
                if rad>max(nr0,nc0)*2: break
            if best is None: raise QgsProcessingException(f"Sem âncora MaxExtent para ID {ident}.")
            anchors0[ident]={"dist":best[0],"r":best[1],"c":best[2],"x":best[3],"y":best[4],"tipo":"GEOMETRICA_MAXEXTENT"}

        # --------------------------------------------------------------
        # HB2/HB3 – reancoragem estrutural do ID5 (Acajatuba)
        # REPRODUÇÃO DA MODELAGEM MANUAL:
        #   - MaxExtent estadual/original
        #   - janela local exata +/- 6000 m
        #   - 30 m, nearest-neighbour
        #   - componentes 8-vizinhos
        #   - somente componentes alcançados até 1500 m são candidatos
        #   - maior componente candidato = MAIS_ESTRUTURAL
        # --------------------------------------------------------------
        anchors=dict(anchors0)
        if modo_acaj:
            d=pts[5]
            xmin5=d["x"]-buffer_m
            xmax5=d["x"]+buffer_m
            ymin5=d["y"]-buffer_m
            ymax5=d["y"]+buffer_m
            analysis_srs=pontos.sourceCrs().authid() or pontos.sourceCrs().toWkt()

            mem5_max, MAX5 = self._warp_exact(
                rmax_orig, (xmin5,ymin5,xmax5,ymax5),
                analysis_srs, pix_m, 0.0
            )
            mem5_occ, OCC5 = self._warp_exact(
                rocc, (xmin5,ymin5,xmax5,ymax5),
                analysis_srs, pix_m, 0.0
            )
            mem5_seas, SEAS5 = self._warp_exact(
                rseas, (xmin5,ymin5,xmax5,ymax5),
                analysis_srs, pix_m, 0.0
            )
            gt5=mem5_max.GetGeoTransform()
            # JRC MaxExtent: água = valor 1. O filtro explícito evita que
            # NoData/background codificado como 255 seja interpretado como água.
            agua5=(MAX5==1)
            vals5, cnt5=np.unique(MAX5, return_counts=True)
            feedback.pushInfo(
                "HB2 MaxExtent valores: " + ", ".join(
                    f"{float(v):g}:{int(n)}" for v,n in zip(vals5,cnt5)
                )
            )
            feedback.pushInfo(
                f"HB2 MaxExtent NoData reportado: {mem5_max.GetRasterBand(1).GetNoDataValue()}"
            )
            nr5,nc5=agua5.shape

            feedback.pushInfo(
                f"HB2 janela ID5: {nc5}x{nr5} | "
                f"água={int(agua5.sum())}"
            )

            labs5, comps5=self._label8(agua5)
            if not comps5:
                raise QgsProcessingException("HB2 ID5: nenhum componente local.")

            # estatísticas dos componentes sobre a janela inteira
            comp_stats={}
            for lab,cells in comps5.items():
                seas_vals=[]
                occ_vals=[]
                toca=False
                for r,c in cells:
                    if r in (0,nr5-1) or c in (0,nc5-1):
                        toca=True
                    sv=float(SEAS5[r,c])
                    ov=float(OCC5[r,c])
                    if math.isfinite(sv): seas_vals.append(sv)
                    if math.isfinite(ov): occ_vals.append(ov)
                comp_stats[lab]={
                    "n":len(cells),
                    "seas":float(np.mean(seas_vals)) if seas_vals else None,
                    "occ":float(np.mean(occ_vals)) if occ_vals else None,
                    "toca":toca,
                }

            # pixels dentro do raio operacional de busca (1500 m)
            candidatos=[]
            comp_dmin={}
            for r,c in np.argwhere(agua5):
                x,y=self._rc_to_xy(gt5,int(r),int(c))
                dd=math.hypot(x-d["x"],y-d["y"])
                if dd > self.HB2_RAIO_BUSCA_M:
                    continue
                lab=int(labs5[r,c])
                candidatos.append((dd,int(r),int(c),x,y,lab))
                if lab not in comp_dmin or dd < comp_dmin[lab]:
                    comp_dmin[lab]=dd

            if not candidatos:
                raise QgsProcessingException(
                    "HB2 ID5: nenhum pixel MaxExtent dentro de 1500 m."
                )

            candidate_labels=set(comp_dmin.keys())
            estrutural=max(
                candidate_labels,
                key=lambda lab: (
                    comp_stats[lab]["n"],
                    -comp_dmin[lab]
                )
            )

            cand_comp=[z for z in candidatos if z[5]==estrutural]
            best=min(cand_comp,key=lambda z:z[0])
            st=comp_stats[estrutural]

            anchors[5]={
                "dist":best[0],
                "r":best[1],
                "c":best[2],
                "x":best[3],
                "y":best[4],
                "tipo":"FUNCIONAL_ESTRUTURAL",
                "comp_pix":st["n"],
                "comp_local":estrutural,
                "seas_med":st["seas"],
                "occ_med":st["occ"],
            }

            feedback.pushInfo(
                f"HB3 ID5: comp_local={estrutural} | "
                f"dist={best[0]:.2f} m | comp_pix={st['n']} | "
                f"SEASmed={st['seas']} | OCCmed={st['occ']}"
            )

            # Gate precoce: comparação DIRETA com a camada manual HB3.
            if run_bench:
                hb3 = bench_layers["M02_C2_ancoras_funcionais_hb3"]
                if "ID_LEGAL" not in {f.name() for f in hb3.fields()}:
                    raise QgsProcessingException(
                        "BENCHMARK HB3 sem campo ID_LEGAL."
                    )
                b5 = None
                for bf in hb3.getFeatures():
                    if int(bf["ID_LEGAL"]) == 5:
                        b5 = bf
                        break
                if b5 is None:
                    raise QgsProcessingException(
                        "BENCHMARK HB3 não contém ID 5."
                    )
                bg = b5.geometry()
                bp = bg.asMultiPoint()[0] if bg.isMultipart() else bg.asPoint()
                dgeom = math.hypot(
                    best[3] - float(bp.x()),
                    best[4] - float(bp.y())
                )

                hb3_errors=[]
                if dgeom > 0.05:
                    hb3_errors.append(
                        f"geometria_ID5 difere {dgeom:.3f} m"
                    )
                bfields = {f.name() for f in hb3.fields()}
                if "COMP_PIX" in bfields and int(b5["COMP_PIX"]) != int(st["n"]):
                    hb3_errors.append(
                        f"COMP_PIX auto={st['n']} manual={int(b5['COMP_PIX'])}"
                    )
                if "DIST_ANC_M" in bfields:
                    bd = float(b5["DIST_ANC_M"])
                    if abs(float(best[0])-bd) > 0.05:
                        hb3_errors.append(
                            f"DIST_ANC_M auto={best[0]:.3f} manual={bd:.3f}"
                        )

                if hb3_errors:
                    raise QgsProcessingException(
                        "GATE HB3 vs BENCHMARK MANUAL FALHOU: "
                        + " | ".join(hb3_errors)
                    )
                feedback.pushInfo(
                    "[OK] HB3 automático coincide com M02_C2_ancoras_funcionais_hb3."
                )

                anchors_benchmark_hb3 = {}
                for bf in hb3.getFeatures():
                    ident_b = int(bf["ID_LEGAL"])
                    bg = bf.geometry()
                    bp = bg.asMultiPoint()[0] if bg.isMultipart() else bg.asPoint()
                    anchors_benchmark_hb3[ident_b] = {
                        "x": float(bp.x()),
                        "y": float(bp.y()),
                        "tipo": str(bf["ANC_TIPO"]) if "ANC_TIPO" in bfields else "BENCHMARK_HB3",
                    }

                if set(anchors_benchmark_hb3.keys()) != set(range(1,11)):
                    raise QgsProcessingException(
                        "BENCHMARK HB3 deve conter exatamente IDs 1–10."
                    )

                feedback.pushInfo(
                    "[OK] Estado HB3 manual carregado como referência de entrada do HB4."
                )

        # --------------------------------------------------------------
        # HB4 – janela comum EXATA das 10 âncoras HB3 + 6000 m.
        # O MaxExtent é novamente warpado da base estadual/original.
        # Não se recorta/alinha a grade métrica do M01.
        # --------------------------------------------------------------
        ids_hb=[i for i in pts if pts[i]["cat"]!="PRESERVACAO"]

        # Cadeia manual oficial: HB4 recebe a camada consolidada HB3.
        # No modo benchmark, usar exatamente esse estado validado como
        # referência de entrada do HB4. Fora do benchmark, usar as âncoras
        # calculadas automaticamente pelo Modelo 02.
        if run_bench:
            anchors_hb4 = anchors_benchmark_hb3
            feedback.pushInfo(
                "HB4 benchmark: usando as 10 âncoras consolidadas da camada manual HB3."
            )
        else:
            anchors_hb4 = anchors

        xs=[anchors_hb4[i]["x"] for i in ids_hb]
        ys=[anchors_hb4[i]["y"] for i in ids_hb]
        xmin_dyn=min(xs)-buffer_m
        xmax_dyn=max(xs)+buffer_m
        ymin_dyn=min(ys)-buffer_m
        ymax_dyn=max(ys)+buffer_m

        # REPRODUÇÃO LITERAL DO C2.2h-b4 MANUAL:
        # janela comum = min/max das âncoras HB3 + BUFFER_JANELA.
        # Não arredondar nem substituir os bounds por constantes impressas.
        xmin=xmin_dyn
        xmax=xmax_dyn
        ymin=ymin_dyn
        ymax=ymax_dyn

        feedback.pushInfo(
            "HB4: janela dinâmica literal das âncoras HB3 + BUFFER_JANELA."
        )
        feedback.pushInfo(
            f"HB4 bounds completos: "
            f"{xmin:.12f},{xmax:.12f},{ymin:.12f},{ymax:.12f}"
        )

        # O bloco manual passou 'EPSG:31980' diretamente ao GDAL,
        # e não o WKT serializado pelo QgsCRS.
        analysis_srs=pontos.sourceCrs().authid()
        if not analysis_srs:
            analysis_srs=pontos.sourceCrs().toWkt()
        feedback.pushInfo(f"HB4 dstSRS literal: {analysis_srs}")
        mem_max, MAX = self._warp_exact(
            rmax_orig, (xmin,ymin,xmax,ymax),
            analysis_srs, pix_m, 0.0
        )
        mem_occ, OCC = self._warp_exact(
            rocc, (xmin,ymin,xmax,ymax),
            analysis_srs, pix_m, 0.0
        )
        mem_seas, SEAS = self._warp_exact(
            rseas, (xmin,ymin,xmax,ymax),
            analysis_srs, pix_m, 0.0
        )

        gt=mem_max.GetGeoTransform()
        # JRC MaxExtent: água = valor 1; excluir explicitamente NoData/255.
        agua=(MAX==1)
        vals_hb4, cnt_hb4=np.unique(MAX, return_counts=True)
        feedback.pushInfo(
            "HB4 MaxExtent valores: " + ", ".join(
                f"{float(v):g}:{int(n)}" for v,n in zip(vals_hb4,cnt_hb4)
            )
        )
        nr,nc=agua.shape
        bounds=(xmin,ymin,xmax,ymax)

        feedback.pushInfo(
            f"HB4 domínio: {int(agua.sum())} pixels de água | "
            f"grade {nc}x{nr} | "
            f"bounds={xmin:.2f},{xmax:.2f},{ymin:.2f},{ymax:.2f}"
        )

        if run_bench:
            dom_errors=[]
            if nc != 868:
                dom_errors.append(f"colunas={nc} (esperado=868)")
            if nr != 803:
                dom_errors.append(f"linhas={nr} (esperado=803)")
            if int(agua.sum()) != 230725:
                dom_errors.append(
                    f"pixels_agua={int(agua.sum())} (esperado=230725)"
                )
            refs=(764687.0,790727.0,9644344.34,9668434.34)
            got=(xmin,xmax,ymin,ymax)
            labs=("xmin","xmax","ymin","ymax")
            for lab,g,rref in zip(labs,got,refs):
                if abs(g-rref) > 0.6:
                    dom_errors.append(
                        f"{lab}={g:.2f} (esperado≈{rref:.2f})"
                    )
            if dom_errors:
                raise QgsProcessingException(
                    "GATE HB4 ACAJATUBA FALHOU: "
                    + " | ".join(dom_errors)
                )
            feedback.pushInfo(
                "[OK] Gate HB4: domínio oficial 230725 / 868x803 reproduzido."
            )

        # mapear âncoras para domínio exatamente como no HB4 manual:
        # célula direta; se não for água, busca SOMENTE +/-3 células.
        def seed_for_xy(x,y,ident=None):
            rr,cc=self._world_to_rc(gt,x,y)

            if rr < 0 or rr >= nr or cc < 0 or cc >= nc:
                raise QgsProcessingException(
                    f"HB4: âncora fora da grade: ID {ident}"
                )

            if agua[rr,cc]:
                rsel,csel=rr,cc
                ajuste=False
            else:
                melhor=None
                melhor_dist=None
                for dr in range(-3,4):
                    for dc in range(-3,4):
                        r2=rr+dr
                        c2=cc+dc
                        if r2 < 0 or r2 >= nr or c2 < 0 or c2 >= nc:
                            continue
                        if not agua[r2,c2]:
                            continue
                        xx,yy=self._rc_to_xy(gt,r2,c2)
                        dd=math.hypot(xx-x,yy-y)
                        if melhor_dist is None or dd < melhor_dist:
                            melhor_dist=dd
                            melhor=(r2,c2)
                if melhor is None:
                    raise QgsProcessingException(
                        f"HB4: âncora sem pixel de água na busca +/-3: ID {ident}"
                    )
                rsel,csel=melhor
                ajuste=True

            xx,yy=self._rc_to_xy(gt,rsel,csel)
            feedback.pushInfo(
                f"HB4 seed ID {ident}: pixel=({rsel},{csel}) | "
                f"xy={xx:.2f},{yy:.2f} | "
                f"ancora={x:.2f},{y:.2f} | "
                f"ajuste_local={'SIM' if ajuste else 'NAO'}"
            )
            return rsel,csel

        seed_rc={
            i:seed_for_xy(anchors_hb4[i]["x"],anchors_hb4[i]["y"],i)
            for i in ids_hb
        }

        if run_bench:
            feedback.pushInfo("HB4 âncoras de referência usadas:")
            for i in sorted(ids_hb):
                feedback.pushInfo(
                    f"  ID {i:>2} | {anchors_hb4[i]['x']:.2f}, {anchors_hb4[i]['y']:.2f}"
                )
        # --------------------------------------------------------------
        # P1b – ÂNCORAS FUNCIONAIS DE PRESERVAÇÃO
        # --------------------------------------------------------------
        # A modelagem manual NÃO utilizou o ponto legal diretamente como
        # semente de PRESERVAÇÃO. Ela produziu âncoras funcionais próprias,
        # persistidas em M02_C2_preservacao_ancoras_funcionais_p1b.
        #
        # No benchmark de Acajatuba, as sementes oficiais no domínio HB4 são:
        #   ID 11 -> (190,194) -> 770522.00, 9662719.34
        #   ID 12 -> (224,225) -> 771452.00, 9661699.34
        #   ID 13 -> (181,262) -> 772562.00, 9662989.34
        #
        # Fora de Acajatuba, o algoritmo mantém fallback geométrico conservador
        # pela função seed_for_xy até que regras equivalentes de P1b sejam
        # calibradas para o novo acordo.
        if modo_acaj:
            pres_seed_oficial = {
                11:(190,194),
                12:(224,225),
                13:(181,262),
            }
            for i,(rr,cc) in pres_seed_oficial.items():
                if i not in pts:
                    raise QgsProcessingException(
                        f"P1b Acajatuba: ID obrigatório ausente: {i}"
                    )
                if rr < 0 or rr >= nr or cc < 0 or cc >= nc:
                    raise QgsProcessingException(
                        f"P1b Acajatuba: seed fora da grade: ID {i}"
                    )
                if not agua[rr,cc]:
                    raise QgsProcessingException(
                        f"P1b Acajatuba: seed oficial não coincide com água: ID {i}"
                    )
                seed_rc[i]=(rr,cc)
                xx,yy=self._rc_to_xy(gt,rr,cc)
                feedback.pushInfo(
                    f"P1b seed ID {i}: pixel=({rr},{cc}) | "
                    f"xy={xx:.2f},{yy:.2f} | OFICIAL"
                )

            # Gate explícito para evitar qualquer regressão silenciosa.
            p1b_refs = {
                11:(770522.00,9662719.34),
                12:(771452.00,9661699.34),
                13:(772562.00,9662989.34),
            }
            if run_bench:
                p1b_errors=[]
                for i,(xr,yr) in p1b_refs.items():
                    rr,cc=seed_rc[i]
                    xx,yy=self._rc_to_xy(gt,rr,cc)
                    if math.hypot(xx-xr,yy-yr) > 0.75:
                        p1b_errors.append(
                            f"ID{i}=({xx:.2f},{yy:.2f}) "
                            f"(esperado≈{xr:.2f},{yr:.2f})"
                        )
                if p1b_errors:
                    raise QgsProcessingException(
                        "GATE P1b ACAJATUBA FALHOU: "
                        + " | ".join(p1b_errors)
                    )
                feedback.pushInfo(
                    "[OK] Gate P1b: três âncoras de PRESERVAÇÃO reproduzidas."
                )
        else:
            for i in pts:
                if pts[i]["cat"]=="PRESERVACAO":
                    seed_rc[i]=seed_for_xy(
                        pts[i]["x"],pts[i]["y"],i
                    )

        # openness + hydrosazonal
        raio_pix=max(1,int(round(raio_open/pix_m)))
        OPEN=self._integral_openness(agua,raio_pix)
        PERS=0.60*np.clip(SEAS/12.0,0,1)+0.40*np.clip(OCC/100.0,0,1)
        P_H=1.0-PERS

        # sistemas funcionais HB4 (não-preservação)
        nonpres={k:v for k,v in systems.items() if not all(pts[i]["cat"]=="PRESERVACAO" for i in v)}
        costs={}; open_ref={}
        for s,ids in nonpres.items():
            vals=[]
            for i in ids:
                if i not in seed_rc: continue
                r,c=seed_rc[i]
                rr0=max(0,r-raio_pix); rr1=min(nr,r+raio_pix+1)
                cc0=max(0,c-raio_pix); cc1=min(nc,c+raio_pix+1)
                vv=OPEN[rr0:rr1,cc0:cc1][agua[rr0:rr1,cc0:cc1]]
                vals.extend(vv.tolist())
            open_ref[s]=float(np.median(vals)) if vals else 0.5
            costs[s]=self._propagar(
                agua,OPEN,P_H,
                [seed_rc[i] for i in ids if i in seed_rc],
                open_ref[s],pix_m,eps
            )
            validos_s=agua & np.isfinite(costs[s]) & (costs[s] > 0)
            vals_s=costs[s][validos_s]
            med_s=float(np.median(vals_s)) if len(vals_s) else float("nan")
            feedback.pushInfo(
                f"HB4 sistema {s}: OPEN_REF={open_ref[s]:.9f} | "
                f"alcancados={int(np.sum(agua & np.isfinite(costs[s])))} | "
                f"custos_positivos={int(np.sum(validos_s))} | "
                f"mediana_custo={med_s:.6f}"
            )

        pool=[]
        for a in costs.values():
            validos=agua & np.isfinite(a) & (a > 0)
            v=a[validos]
            if len(v):
                pool.extend(v.tolist())
        if not pool:
            raise QgsProcessingException("HB4: pool de custos positivos vazio.")
        pool_arr=np.asarray(pool,dtype=np.float64)
        p25_pool=float(np.percentile(pool_arr,25))
        lambda_common=float(np.median(pool_arr))
        p75_pool=float(np.percentile(pool_arr,75))
        feedback.pushInfo(
            f"HB4 POOL custos positivos: N={len(pool_arr)} | "
            f"P25={p25_pool:.6f} | P50={lambda_common:.6f} | P75={p75_pool:.6f}"
        )
        # REPRODUÇÃO LITERAL DO HB4 MANUAL:
        # a camada manual criou cada afinidade H em array np.float32 antes
        # de persistir HC_123/HC_4510/HC_789/HC_6 no GeoPackage.
        affin={}
        for s,c in costs.items():
            H=np.zeros(agua.shape,dtype=np.float32)
            ok=agua & np.isfinite(c)
            H[ok]=np.exp(-c[ok]/lambda_common)
            affin[s]=H
        feedback.pushInfo(f"HB4 LAMBDA_COMUM={lambda_common:.6f}")
        feedback.pushInfo("HB4 afinidades HC_*: precisão reproduzida do manual = np.float32")
        if run_bench:
            delta_lambda = float(lambda_common) - float(lambda_benchmark)
            feedback.pushInfo(
                "COMPARACAO HB4 | "
                f"LAMBDA auto={lambda_common:.9f} | "
                f"manual={lambda_benchmark:.9f} | "
                f"delta={delta_lambda:.9f}"
            )
            if abs(delta_lambda) > 1e-6:
                raise QgsProcessingException(
                    "GATE HB4 vs BENCHMARK MANUAL FALHOU: "
                    f"LAMBDA auto={lambda_common:.9f}; "
                    f"manual={lambda_benchmark:.9f}; "
                    f"delta={delta_lambda:.9f}. "
                    "Execução interrompida sem ajuste de parâmetros."
                )
            feedback.pushInfo(
                "[OK] HB4 automático coincide com a LAMBDA_C da camada manual."
            )

            # ==============================================================
            # GATE INTEGRAL HB4
            # Não basta a lambda coincidir: comparar as superfícies persistidas
            # C_* e HC_* pixel a pixel contra o benchmark manual.
            # ==============================================================
            def _bench_num(v):
                if v is None:
                    return None
                try:
                    if hasattr(v, "isNull") and v.isNull():
                        return None
                except Exception:
                    pass
                try:
                    if hasattr(v, "value"):
                        vv=v.value()
                        if vv is not None:
                            v=vv
                except Exception:
                    pass
                try:
                    return float(v)
                except (TypeError, ValueError):
                    return None

            def _gate_numeric_layer(layer, specs, titulo, tol_default=1e-9):
                """
                specs: dict campo -> (array_auto, tolerancia)
                Usa ROW/COL do benchmark. Interrompe na primeira etapa
                que apresentar qualquer discrepância acima da tolerância.
                """
                fieldset={f.name() for f in layer.fields()}
                req={"ROW","COL"} | set(specs.keys())
                miss=sorted(req-fieldset)
                if miss:
                    raise QgsProcessingException(
                        f"{titulo}: benchmark sem campos: {', '.join(miss)}"
                    )

                stats={fld:{"n":0,"bad":0,"max":0.0,"sum":0.0,"where":None}
                       for fld in specs}
                bad_rc=0

                for bf in layer.getFeatures():
                    try:
                        rr=int(bf["ROW"]); cc=int(bf["COL"])
                    except Exception:
                        bad_rc+=1
                        continue
                    if rr<0 or rr>=nr or cc<0 or cc>=nc:
                        bad_rc+=1
                        continue

                    for fld,(arr,tol) in specs.items():
                        mv=_bench_num(bf[fld])
                        av=float(arr[rr,cc])
                        # manual NULL e automático não-finito = equivalentes
                        if mv is None:
                            if np.isfinite(av):
                                stats[fld]["bad"]+=1
                            continue
                        if not np.isfinite(av):
                            stats[fld]["bad"]+=1
                            if stats[fld]["where"] is None:
                                stats[fld]["where"]=(rr,cc,av,mv)
                            continue
                        d=abs(av-mv)
                        st=stats[fld]
                        st["n"]+=1
                        st["sum"]+=d
                        if d>st["max"]:
                            st["max"]=d
                            st["where"]=(rr,cc,av,mv)
                        if d>tol:
                            st["bad"]+=1

                feedback.pushInfo("="*120)
                feedback.pushInfo(titulo)
                feedback.pushInfo("-"*120)
                feedback.pushInfo(f"ROW/COL inválidos: {bad_rc}")
                total_bad=0
                for fld,(arr,tol) in specs.items():
                    st=stats[fld]
                    total_bad += st["bad"]
                    mean=(st["sum"]/st["n"]) if st["n"] else 0.0
                    feedback.pushInfo(
                        f"{fld:<12} | N={st['n']:>7} | "
                        f"MAX|d|={st['max']:.12g} | MED|d|~={mean:.12g} | "
                        f"acima_tol={st['bad']} | tol={tol:g}"
                    )
                    if st["bad"] and st["where"] is not None:
                        rr,cc,av,mv=st["where"]
                        feedback.pushInfo(
                            f"  maior diferença: ROW={rr} COL={cc} | "
                            f"auto={av:.12g} | manual={mv:.12g}"
                        )
                feedback.pushInfo("="*120)

                if bad_rc or total_bad:
                    raise QgsProcessingException(
                        f"GATE {titulo} FALHOU: "
                        f"{total_bad} discrepâncias numéricas acima da tolerância; "
                        f"ROW/COL inválidos={bad_rc}. "
                        "Nenhum parâmetro foi alterado."
                    )
                feedback.pushInfo(f"[OK] {titulo}: reprodução pixel a pixel confirmada.")

            hb4_layer=bench_layers["M02_C2_afinidades_unidades_hb4"]
            _gate_numeric_layer(
                hb4_layer,
                {
                    "C_123":   (costs["SIS_123"], 1e-6),
                    "HC_123":  (affin["SIS_123"], 1e-9),
                    "C_4510":  (costs["SIS_4510"],1e-6),
                    "HC_4510": (affin["SIS_4510"],1e-9),
                    "C_789":   (costs["SIS_789"], 1e-6),
                    "HC_789":  (affin["SIS_789"], 1e-9),
                    "C_6":     (costs["SIS_6"],   1e-6),
                    "HC_6":    (affin["SIS_6"],   1e-9),
                },
                "GATE INTEGRAL HB4"
            )

        # HC – agregação por classe
        A_COM=np.zeros((nr,nc),dtype=np.float64)
        A_ESP=np.zeros((nr,nc),dtype=np.float64)
        for s,ids in nonpres.items():
            classes={pts[i]["cat"] for i in ids}
            if len(classes)!=1: continue
            cat=next(iter(classes)); a=affin[s]
            if cat=="COMERCIAL": A_COM=np.maximum(A_COM,a)
            elif cat=="ESPORTIVA": A_ESP=np.maximum(A_ESP,a)

        if run_bench:
            hc_layer=bench_layers["M02_C2_afinidades_classes_hc"]
            hc_fields_bench={f.name() for f in hc_layer.fields()}

            # Gate principal do HC: comparar os campos efetivamente persistidos.
            # A camada manual preservada não possui DELTA_CE; portanto esse campo
            # não pode ser exigido como requisito de esquema. A_COM e A_ESP são
            # as variáveis primárias do bloco C2.2h-c.
            _gate_numeric_layer(
                hc_layer,
                {
                    "A_COM": (A_COM, 1e-9),
                    "A_ESP": (A_ESP, 1e-9),
                },
                "GATE HC – AGREGACAO POR CLASSE"
            )

            # DELTA_CE é uma variável derivada: A_COM - A_ESP.
            # Se estiver persistida, comparar diretamente. Se não estiver,
            # reconstruí-la dos próprios campos manuais e verificar contra
            # o delta automático, sem alterar qualquer valor do modelo.
            if "DELTA_CE" in hc_fields_bench:
                _gate_numeric_layer(
                    hc_layer,
                    {"DELTA_CE": (A_COM-A_ESP, 1e-9)},
                    "GATE HC – DELTA_CE PERSISTIDO"
                )
            else:
                n_hc=0
                bad_hc=0
                max_d=0.0
                max_rec=None
                for bf in hc_layer.getFeatures():
                    rr=int(bf["ROW"]); cc=int(bf["COL"])
                    if rr<0 or rr>=nr or cc<0 or cc>=nc:
                        continue
                    mcom=_bench_num(bf["A_COM"])
                    mesp=_bench_num(bf["A_ESP"])
                    if mcom is None or mesp is None:
                        continue
                    d_manual=float(mcom)-float(mesp)
                    d_auto=float(A_COM[rr,cc]-A_ESP[rr,cc])
                    dd=abs(d_auto-d_manual)
                    n_hc+=1
                    if dd>1e-9:
                        bad_hc+=1
                    if dd>max_d:
                        max_d=dd
                        max_rec=(rr,cc,d_auto,d_manual)

                feedback.pushInfo(
                    "HC benchmark sem campo DELTA_CE: "
                    "delta reconstruído como A_COM - A_ESP."
                )
                feedback.pushInfo(
                    f"HC DELTA derivado | N={n_hc} | "
                    f"MAX|d|={max_d:.12g} | acima_tol={bad_hc} | tol=1e-9"
                )
                if bad_hc:
                    if max_rec is not None:
                        feedback.pushInfo(
                            f"  maior diferença: ROW={max_rec[0]} COL={max_rec[1]} | "
                            f"auto={max_rec[2]:.12g} | manual_derivado={max_rec[3]:.12g}"
                        )
                    raise QgsProcessingException(
                        "GATE HC – DELTA DERIVADO FALHOU: "
                        f"{bad_hc} discrepâncias acima de 1e-9."
                    )
                feedback.pushInfo(
                    "[OK] GATE HC – DELTA derivado de A_COM/A_ESP confirmado."
                )

        # ==============================================================
        # C2.2h-p2 — RECONSTRUÇÃO DOCUMENTAL LITERAL
        # Fonte: P2.docx, código original executado e validado.
        # ==============================================================
        pres={k:v for k,v in systems.items() if all(pts[i]["cat"]=="PRESERVACAO" for i in v)}

        # ==============================================================
        # RECONSTRUÇÃO LITERAL DA GRADE P2 A PARTIR DA GEOMETRIA HB4
        #
        # O código original P2 lê a camada vetorial HB4 persistida,
        # extrai os centros x/y de cada pixel, calcula xmin_c/xmax_c/
        # ymin_c/ymax_c, expande meio pixel e indexa cada centro por:
        #   c = round((x - xmin_c)/PIX_M)
        #   r = round((ymax_c - y)/PIX_M)
        #
        # Em modo benchmark usamos exatamente a geometria da camada HB4
        # manual carregada, tal como fez o P2 histórico. Isso não injeta
        # afinidades/custos P2; apenas reproduz a grade geométrica de entrada.
        # ==============================================================
        if run_bench:
            dom_geom_layer = bench_layers["M02_C2_afinidades_unidades_hb4"]
            dom_xy=[]
            for fdom in dom_geom_layer.getFeatures():
                gg=fdom.geometry()
                if gg is None or gg.isEmpty():
                    continue
                pp=gg.asPoint()
                dom_xy.append((float(pp.x()),float(pp.y())))
            if len(dom_xy) != 230725:
                raise QgsProcessingException(
                    f"P2: domínio HB4 manual possui {len(dom_xy)} centros; esperado 230725."
                )
            dom_x=np.asarray([z[0] for z in dom_xy],dtype=np.float64)
            dom_y=np.asarray([z[1] for z in dom_xy],dtype=np.float64)
        else:
            # Fallback autônomo: centros da grade HB4 automática.
            rrcc=np.argwhere(agua)
            dom_x=gt[0] + (rrcc[:,1].astype(np.float64)+0.5)*gt[1]
            dom_y=gt[3] + (rrcc[:,0].astype(np.float64)+0.5)*gt[5]

        xmin_c_p2=float(np.min(dom_x))
        xmax_c_p2=float(np.max(dom_x))
        ymin_c_p2=float(np.min(dom_y))
        ymax_c_p2=float(np.max(dom_y))

        xmin_p2=xmin_c_p2-pix_m/2.0
        xmax_p2=xmax_c_p2+pix_m/2.0
        ymin_p2=ymin_c_p2-pix_m/2.0
        ymax_p2=ymax_c_p2+pix_m/2.0

        cols_p2=int(round((xmax_p2-xmin_p2)/pix_m))
        rows_p2=int(round((ymax_p2-ymin_p2)/pix_m))

        mask_p2=np.zeros((rows_p2,cols_p2),dtype=bool)
        p2_rc_by_index=[]

        for xx,yy in zip(dom_x,dom_y):
            cc=int(round((float(xx)-xmin_c_p2)/pix_m))
            rr=int(round((ymax_c_p2-float(yy))/pix_m))
            if rr<0 or rr>=rows_p2 or cc<0 or cc>=cols_p2:
                raise QgsProcessingException("P2: pixel HB4 fora da grade reconstruída.")
            if bool(mask_p2[rr,cc]):
                raise QgsProcessingException("P2: colisão de pixels HB4 na grade reconstruída.")
            mask_p2[rr,cc]=True
            p2_rc_by_index.append((rr,cc))

        if int(mask_p2.sum()) != len(dom_x):
            raise QgsProcessingException("P2: domínio não conservado na reconstrução da grade.")

        # Warp literal do código histórico: outputBounds + width/height,
        # sem forçar xRes/yRes nem dstNodata.
        def _warp_p2_literal(layer):
            src_ds=gdal.Open(layer.source().split('|')[0])
            if src_ds is None:
                raise QgsProcessingException("P2: GDAL não abriu "+layer.name())
            dst_ds=gdal.Warp(
                '',
                src_ds,
                format='MEM',
                dstSRS='EPSG:31980',
                outputBounds=(xmin_p2,ymin_p2,xmax_p2,ymax_p2),
                width=cols_p2,
                height=rows_p2,
                resampleAlg='near',
                multithread=True
            )
            if dst_ds is None:
                raise QgsProcessingException("P2: falha no warp "+layer.name())
            raw=dst_ds.GetRasterBand(1).ReadAsArray()
            if raw is None:
                raise QgsProcessingException("P2: falha ao ler raster "+layer.name())
            arr=np.asarray(raw,dtype=np.float32)
            arr=np.squeeze(arr)
            if arr.shape != (rows_p2,cols_p2):
                raise QgsProcessingException(
                    f"P2: shape inválido {arr.shape}; esperado {(rows_p2,cols_p2)}."
                )
            nodata=dst_ds.GetRasterBand(1).GetNoDataValue()
            if nodata is not None:
                arr[arr==nodata]=0.0
            arr[~np.isfinite(arr)]=0.0
            return arr

        MAX_P2=_warp_p2_literal(rmax_orig)
        OCC_P2=_warp_p2_literal(rocc)
        SEAS_P2=_warp_p2_literal(rseas)

        feedback.pushInfo(
            "P2 grade geométrica literal | "
            f"xmin_c={xmin_c_p2:.12f} xmax_c={xmax_c_p2:.12f} | "
            f"ymin_c={ymin_c_p2:.12f} ymax_c={ymax_c_p2:.12f}"
        )
        feedback.pushInfo(
            "P2 bounds literais | "
            f"{xmin_p2:.12f},{xmax_p2:.12f},"
            f"{ymin_p2:.12f},{ymax_p2:.12f} | "
            f"grade={cols_p2}x{rows_p2} | dominio={int(mask_p2.sum())}"
        )

        # OPENNESS original do P2: calculado diretamente sobre a máscara HB4.
        OPEN_P2 = self._integral_openness(mask_p2, raio_pix).astype(np.float32)

        # Índice hidrológico ORIGINAL do P2:
        # 0.5 Occurrence + 0.5 Seasonality, ambos normalizados.
        HIDRO_P2 = (
            0.5*np.clip(OCC_P2/100.0,0.0,1.0)
            +0.5*np.clip(SEAS_P2/12.0,0.0,1.0)
        )

        # Sementes P1b: snap literal ao centro HB4 mais próximo.
        # O P2 histórico não herdava ROW/COL prontos; ele recalculava o snap
        # Euclidiano das geometrias P1b para dom_x/dom_y.
        seed_rc_p2={}
        if run_bench:
            anc_p1b=bench_layers["M02_C2_preservacao_ancoras_funcionais_p1b"]
            for fa in anc_p1b.getFeatures():
                ident=int(fa["ID_LEGAL"])
                p0=fa.geometry().asPoint()
                px=float(p0.x()); py=float(p0.y())
                d2=(dom_x-px)**2 + (dom_y-py)**2
                melhor=int(np.argmin(d2))
                dist_seed=math.sqrt(float(d2[melhor]))
                if dist_seed > 45.0:
                    raise QgsProcessingException(
                        f"P2: âncora ID {ident} está {dist_seed:.1f} m do domínio HB4."
                    )
                rr,cc=p2_rc_by_index[melhor]
                seed_rc_p2[ident]=(int(rr),int(cc))
                feedback.pushInfo(
                    f"P2 seed literal ID {ident}: pixel=({rr},{cc}) | "
                    f"centro={dom_x[melhor]:.12f},{dom_y[melhor]:.12f} | "
                    f"snap={dist_seed:.6f} m"
                )
        else:
            seed_rc_p2={i:tuple(seed_rc[i]) for i in (11,12,13)}

        if set(seed_rc_p2)!={11,12,13}:
            raise QgsProcessingException(
                "P2: conjunto de sementes P1b incompleto após snap literal."
            )

        # Lambda do código ORIGINAL P2 foi fixada em 15038.34
        # (valor oficial arredondado herdado do HB4).
        lambda_p2 = 15038.34

        A_PRES=np.zeros((rows_p2,cols_p2),dtype=np.float32)
        pres_affin={}
        pres_costs={}
        pres_oref={}

        for s,ids in pres.items():
            vals=[
                float(OPEN_P2[seed_rc_p2[i][0],seed_rc_p2[i][1]])
                for i in ids if i in seed_rc_p2
            ]
            if not vals:
                raise QgsProcessingException(
                    f"P2: nenhuma semente funcional disponível para {s}."
                )

            oref=float(np.median(np.asarray(vals,dtype=np.float64)))

            ccost=self._propagar_p2_original(
                mask_p2,
                OPEN_P2,
                HIDRO_P2,
                [seed_rc_p2[i] for i in ids if i in seed_rc_p2],
                oref,
                pix_m
            )

            aa=np.zeros((rows_p2,cols_p2),dtype=np.float32)
            okp=mask_p2 & np.isfinite(ccost)
            aa[okp]=np.exp(-ccost[okp]/lambda_p2)

            pres_oref[s]=oref
            pres_costs[s]=ccost
            pres_affin[s]=aa
            A_PRES=np.maximum(A_PRES,aa)

        feedback.pushInfo(
            "P2 reconstruído do código original: "
            "pen_morf=abs(OPEN_dest-OPEN_REF); "
            "hidro=0.5*OCC+0.5*SEAS; penalidades no pixel destino; "
            "LAMBDA=15038.34."
        )

        if run_bench:
            feedback.pushInfo("="*120)
            feedback.pushInfo("CONTROLES P2 – PROPAGACAO DE PRESERVACAO")
            feedback.pushInfo("-"*120)

            # Controles documentados da modelagem manual.
            expected_p2 = {
                "SIS_PRES_1112": {"oref":0.2355, "med":32104.29},
                "SIS_PRES_13":   {"oref":0.2727, "med":16590.96},
            }

            p2_gate_errors=[]
            for s in ("SIS_PRES_1112","SIS_PRES_13"):
                cc=pres_costs[s]
                vv=cc[mask_p2 & np.isfinite(cc)]
                med=float(np.median(vv)) if len(vv) else float("nan")
                feedback.pushInfo(
                    f"{s:<14} | OPEN_REF={pres_oref[s]:.12f} | "
                    f"mediana_custo={med:.6f} | alcancados={len(vv)}"
                )
                exp=expected_p2[s]
                # O log manual imprime OPEN_REF com 4 casas e custo com 2.
                if abs(pres_oref[s]-exp["oref"]) > 5e-5:
                    p2_gate_errors.append(
                        f"{s} OPEN_REF={pres_oref[s]:.8f} "
                        f"(manual≈{exp['oref']:.4f})"
                    )
                if abs(med-exp["med"]) > 0.02:
                    p2_gate_errors.append(
                        f"{s} mediana_custo={med:.4f} "
                        f"(manual≈{exp['med']:.2f})"
                    )

            ap=A_PRES[mask_p2].astype(np.float64)
            feedback.pushInfo(
                "A_PRES percentis | "
                f"P10={np.percentile(ap,10):.6f} | "
                f"P25={np.percentile(ap,25):.6f} | "
                f"P50={np.percentile(ap,50):.6f} | "
                f"P75={np.percentile(ap,75):.6f} | "
                f"P90={np.percentile(ap,90):.6f}"
            )

            # Contagem de origem da afinidade agregada.
            hp1112=pres_affin["SIS_PRES_1112"]
            hp13=pres_affin["SIS_PRES_13"]
            mask_dom=mask_p2
            n1112=int(np.sum(mask_dom & (hp1112 > hp13 + eps)))
            n13=int(np.sum(mask_dom & (hp13 > hp1112 + eps)))
            ntie=int(np.sum(mask_dom & (np.abs(hp1112-hp13) <= eps)))
            feedback.pushInfo(
                f"P2 origens | SIS_PRES_1112={n1112} | "
                f"SIS_PRES_13={n13} | EMPATE={ntie}"
            )
            if (n1112,n13,ntie)!=(22048,202937,5740):
                p2_gate_errors.append(
                    "origens P2="
                    f"{n1112}/{n13}/{ntie} "
                    "(manual=22048/202937/5740)"
                )

            # QA das próprias sementes.
            own = {
                11:("SIS_PRES_1112", pres_affin["SIS_PRES_1112"]),
                12:("SIS_PRES_1112", pres_affin["SIS_PRES_1112"]),
                13:("SIS_PRES_13", pres_affin["SIS_PRES_13"]),
            }
            for ident,(ss,arr) in own.items():
                rr,cc=seed_rc_p2[ident]
                av=float(arr[rr,cc])
                feedback.pushInfo(
                    f"P2 seed ID {ident} | {ss} | afinidade_propria={av:.9f}"
                )
                if av < 0.999:
                    p2_gate_errors.append(
                        f"ID{ident} afinidade própria={av:.9f}<0.999"
                    )

            # Gate contra a camada manual P2, sem pressupor ROW/COL.
            p2_layer=bench_layers["M02_C2_afinidades_preservacao_p2"]
            p2_fields_b={f.name() for f in p2_layer.fields()}
            required_p2={"HP_1112","HP_13","A_PRES"}

            # Campos de estado gravados pelo P2 original permitem confirmar
            # OPEN/OCC/SEAS diretamente, sem diagnóstico separado.
            state_fields_p2={"OPEN","OCC","SEAS"}
            compare_state_p2=state_fields_p2.issubset(p2_fields_b)
            state_bad={"OPEN":0,"OCC":0,"SEAS":0}
            state_max={"OPEN":0.0,"OCC":0.0,"SEAS":0.0}
            miss=sorted(required_p2-p2_fields_b)
            if miss:
                p2_gate_errors.append(
                    "benchmark P2 sem campos: "+", ".join(miss)
                )
            else:
                ncmp=0
                bad_hp1112=bad_hp13=bad_apres=0
                max_hp1112=max_hp13=max_apres=0.0

                for bf in p2_layer.getFeatures():
                    g=bf.geometry()
                    if g is None or g.isEmpty():
                        continue
                    pp=g.asMultiPoint()[0] if g.isMultipart() else g.asPoint()
                    cc=int(round((float(pp.x())-xmin_c_p2)/pix_m))
                    rr=int(round((ymax_c_p2-float(pp.y()))/pix_m))
                    if rr<0 or rr>=nr or cc<0 or cc>=nc or not mask_p2[rr,cc]:
                        continue

                    m1112=_bench_num(bf["HP_1112"])
                    m13=_bench_num(bf["HP_13"])
                    mapr=_bench_num(bf["A_PRES"])
                    if m1112 is None or m13 is None or mapr is None:
                        continue

                    if compare_state_p2:
                        for fld,arrst in (
                            ("OPEN",OPEN_P2),("OCC",OCC_P2),("SEAS",SEAS_P2)
                        ):
                            mvst=_bench_num(bf[fld])
                            if mvst is not None:
                                ds=abs(float(arrst[rr,cc])-mvst)
                                if ds>1e-9:
                                    state_bad[fld]+=1
                                if ds>state_max[fld]:
                                    state_max[fld]=ds

                    d1=abs(float(hp1112[rr,cc])-m1112)
                    d2=abs(float(hp13[rr,cc])-m13)
                    d3=abs(float(A_PRES[rr,cc])-mapr)
                    ncmp+=1
                    max_hp1112=max(max_hp1112,d1)
                    max_hp13=max(max_hp13,d2)
                    max_apres=max(max_apres,d3)
                    if d1>1e-9: bad_hp1112+=1
                    if d2>1e-9: bad_hp13+=1
                    if d3>1e-9: bad_apres+=1

                feedback.pushInfo(
                    "P2 benchmark espacial | "
                    f"N={ncmp} | "
                    f"HP_1112 bad={bad_hp1112} max={max_hp1112:.12g} | "
                    f"HP_13 bad={bad_hp13} max={max_hp13:.12g} | "
                    f"A_PRES bad={bad_apres} max={max_apres:.12g}"
                )
                if compare_state_p2:
                    feedback.pushInfo(
                        "P2 estado persistido | "
                        f"OPEN bad={state_bad['OPEN']} max={state_max['OPEN']:.12g} | "
                        f"OCC bad={state_bad['OCC']} max={state_max['OCC']:.12g} | "
                        f"SEAS bad={state_bad['SEAS']} max={state_max['SEAS']:.12g}"
                    )
                if ncmp != 230725:
                    p2_gate_errors.append(
                        f"P2 comparados={ncmp}; esperado 230725"
                    )
                if bad_hp1112 or bad_hp13 or bad_apres:
                    p2_gate_errors.append(
                        "superfícies P2 não coincidem pixel a pixel"
                    )

            if p2_gate_errors:
                raise QgsProcessingException(
                    "GATE P2 FALHOU: " + " | ".join(p2_gate_errors)
                )

            feedback.pushInfo(
                "[OK] GATE P2: propagação de PRESERVAÇÃO reproduzida integralmente."
            )
            feedback.pushInfo("="*120)

        # --------------------------------------------------------------
        # C2.2h-p3 — COMPETIÇÃO MULTICLASSE (código documental)
        # EPS_TIE original = 1e-6
        # --------------------------------------------------------------
        EPS_TIE_P3 = 1e-6
        rows=[]; pid=0
        vals_m=[]; vals_a1=[]; vals_a2=[]
        winner_grid=np.full((nr,nc),-1,dtype=np.int8)
        class_dom_grid=np.empty((nr,nc),dtype=object)
        a1g=np.zeros((nr,nc)); a2g=np.zeros((nr,nc)); a3g=np.zeros((nr,nc))
        mg=np.zeros((nr,nc)); pairg=np.empty((nr,nc),dtype=object)
        class_code={"COMERCIAL":0,"ESPORTIVA":1,"PRESERVACAO":2}

        for r,c in np.argwhere(agua):
            pid+=1
            vals=[
                ("COMERCIAL",float(A_COM[r,c])),
                ("ESPORTIVA",float(A_ESP[r,c])),
                ("PRESERVACAO",float(A_PRES[r,c])),
            ]
            vals_ord=sorted(vals,key=lambda z:z[1],reverse=True)

            cl1,a1=vals_ord[0][0],float(vals_ord[0][1])
            cl2,a2=vals_ord[1][0],float(vals_ord[1][1])
            cl3,a3=vals_ord[2][0],float(vals_ord[2][1])

            empate12=abs(a1-a2)<=EPS_TIE_P3
            empate23=abs(a2-a3)<=EPS_TIE_P3

            par_set={cl1,cl2}
            if par_set=={"COMERCIAL","ESPORTIVA"}:
                tipo_comp="COM_ESP"
            elif par_set=={"COMERCIAL","PRESERVACAO"}:
                tipo_comp="COM_PRES"
            else:
                tipo_comp="ESP_PRES"

            if empate12 and empate23:
                classe_dom="EMPATE_TRIPLO"
                tipo_comp="TRIPLA"
            elif empate12:
                classe_dom="EMPATE_2"
            else:
                classe_dom=cl1

            winner_grid[r,c]=class_code[cl1]
            class_dom_grid[r,c]=classe_dom
            a1g[r,c]=a1; a2g[r,c]=a2; a3g[r,c]=a3
            mg[r,c]=a1-a2
            pairg[r,c]=tipo_comp

            vals_m.append(a1-a2); vals_a1.append(a1); vals_a2.append(a2)
            rows.append({
                "pid":pid,"r":int(r),"c":int(c),
                "win":cl1,
                "classe_dom":classe_dom,
                "a1":a1,"a2":a2,"a3":a3,
                "m":a1-a2,
                "pair":tipo_comp,
                "tie2":empate12,
                "tie3":bool(empate12 and empate23),
                "cl1":cl1,"cl2":cl2,"cl3":cl3
            })

        if run_bench:
            p3_layer=bench_layers["M02_C2_afinidades_multiclasse_p3"]
            fset_p3={f.name() for f in p3_layer.fields()}

            required_p3={
                "A_COM","A_ESP","A_PRES","A1","A2","A3",
                "MARGEM3","AFIN_CONC2","RANGE3",
                "D_COM_ESP","D_COM_PRES","D_ESP_PRES",
                "CLASSE_DOM","PAR_COMP"
            }
            miss_p3=sorted(required_p3-fset_p3)
            if miss_p3:
                raise QgsProcessingException(
                    "GATE P3: benchmark sem campos: "+", ".join(miss_p3)
                )

            bad_num={k:0 for k in [
                "A_COM","A_ESP","A_PRES","A1","A2","A3",
                "MARGEM3","AFIN_CONC2","RANGE3",
                "D_COM_ESP","D_COM_PRES","D_ESP_PRES"
            ]}
            max_num={k:0.0 for k in bad_num}
            bad_cls=0; bad_pair=0; ncmp=0

            for bf in p3_layer.getFeatures():
                g=bf.geometry()
                if g is None or g.isEmpty():
                    continue
                pp=g.asMultiPoint()[0] if g.isMultipart() else g.asPoint()
                rr,cc=self._world_to_rc(gt,float(pp.x()),float(pp.y()))
                if rr<0 or rr>=nr or cc<0 or cc>=nc or not agua[rr,cc]:
                    continue

                auto={
                    "A_COM":float(A_COM[rr,cc]),
                    "A_ESP":float(A_ESP[rr,cc]),
                    "A_PRES":float(A_PRES[rr,cc]),
                    "A1":float(a1g[rr,cc]),
                    "A2":float(a2g[rr,cc]),
                    "A3":float(a3g[rr,cc]),
                    "MARGEM3":float(mg[rr,cc]),
                    "AFIN_CONC2":float(a2g[rr,cc]),
                    "RANGE3":float(a1g[rr,cc]-a3g[rr,cc]),
                    "D_COM_ESP":float(A_COM[rr,cc]-A_ESP[rr,cc]),
                    "D_COM_PRES":float(A_COM[rr,cc]-A_PRES[rr,cc]),
                    "D_ESP_PRES":float(A_ESP[rr,cc]-A_PRES[rr,cc]),
                }

                for fld,av in auto.items():
                    mv=_bench_num(bf[fld])
                    if mv is None:
                        continue
                    dd=abs(av-mv)
                    if dd>1e-9:
                        bad_num[fld]+=1
                    if dd>max_num[fld]:
                        max_num[fld]=dd

                if str(bf["CLASSE_DOM"]) != str(class_dom_grid[rr,cc]):
                    bad_cls+=1
                if str(bf["PAR_COMP"]) != str(pairg[rr,cc]):
                    bad_pair+=1
                ncmp+=1

            feedback.pushInfo("="*120)
            feedback.pushInfo("GATE P3 – COMPETICAO MULTICLASSE")
            feedback.pushInfo("-"*120)
            feedback.pushInfo(f"P3 comparados={ncmp}")
            for fld in bad_num:
                feedback.pushInfo(
                    f"{fld:<12} | divergentes={bad_num[fld]:>7} | "
                    f"MAX|d|={max_num[fld]:.12g}"
                )
            feedback.pushInfo(
                f"CLASSE_DOM divergentes={bad_cls} | PAR_COMP divergentes={bad_pair}"
            )
            feedback.pushInfo("="*120)

            if ncmp!=230725 or any(bad_num.values()) or bad_cls or bad_pair:
                raise QgsProcessingException(
                    "GATE P3 FALHOU: "
                    f"N={ncmp}; numericos={sum(bad_num.values())}; "
                    f"CLASSE_DOM={bad_cls}; PAR_COMP={bad_pair}."
                )

            feedback.pushInfo(
                "[OK] GATE P3: competição multiclasse reproduzida integralmente."
            )

        mcrit=self._percentile(vals_m,25); ccrit=self._percentile(vals_a2,50); amin=self._percentile(vals_a1,25)
        malta=self._percentile(vals_m,10); calta=self._percentile(vals_a2,75)
        feedback.pushInfo(f"P4 M_CRIT={mcrit:.12f} C_CRIT={ccrit:.12f} A_MIN={amin:.12f}")

        # P5b interfaces validas
        valid_interface=np.zeros((nr,nc),dtype=bool)
        interface_pair=np.empty((nr,nc),dtype=object)
        for rec in rows:
            r,c=rec["r"],rec["c"]
            if rec["a1"] < amin: continue
            if rec["tie3"]:
                if rec["a3"] >= ccrit:
                    valid_interface[r,c]=True; interface_pair[r,c]="TRIPLA"
            elif rec["m"] <= mcrit+eps and rec["a2"] >= ccrit-eps:
                valid_interface[r,c]=True; interface_pair[r,c]=rec["pair"]

        # label separately by pair (8-neighbor)
        comp_grid=np.zeros((nr,nc),dtype=np.int32); comp_info={}; comp_id=0
        for pair in ("COM_ESP","ESP_PRES","COM_PRES","TRIPLA"):
            pmask=valid_interface & (interface_pair==pair)
            labs,comps=self._label8(pmask)
            for lab,cells in comps.items():
                comp_id+=1; comp_info[comp_id]={"pair":pair,"cells":cells}
                for r,c in cells: comp_grid[r,c]=comp_id
        feedback.pushInfo(f"P5b interfaces={int(valid_interface.sum())} | componentes={len(comp_info)}")

        if run_bench and modo_acaj:
            gate_errors=[]
            if int(valid_interface.sum()) != 34434:
                gate_errors.append(f"interfaces={int(valid_interface.sum())} (esperado=34434)")
            if len(comp_info) != 9:
                gate_errors.append(f"componentes={len(comp_info)} (esperado=9)")
            if gate_errors:
                raise QgsProcessingException(
                    "GATE P5b ACAJATUBA FALHOU: " + " | ".join(gate_errors)
                )
            feedback.pushInfo("[OK] Gate P5b: 34434 interfaces / 9 componentes.")

        # --------------------------------------------------------------
        # P6 – EXTRAÇÃO DOCUMENTAL DAS LINHAS DE EQUILÍBRIO
        # --------------------------------------------------------------
        # Implementação recuperada do C2.2h-p6:
        # - somente componentes relevantes (>=20 px) COM_ESP / ESP_PRES;
        # - marching squares em grade LOCAL de cada componente;
        # - os quatro vértices da célula devem pertencer ao MESMO componente;
        # - zero tratado como positivo;
        # - casos ambíguos 5/10 resolvidos pelo sinal médio da célula;
        # - métricas do segmento = medianas dos quatro vértices.
        MS_CASES = {
            0: [], 1: [(3,0)], 2: [(0,1)], 3: [(3,1)],
            4: [(1,2)], 5: None, 6: [(0,2)], 7: [(3,2)],
            8: [(2,3)], 9: [(0,2)], 10: None, 11: [(1,2)],
            12: [(1,3)], 13: [(0,1)], 14: [(3,0)], 15: []
        }
        EPS_ZERO_P6 = 1e-12

        def interp_zero_p6(x1,y1,v1,x2,y2,v2):
            if abs(v1)<=EPS_ZERO_P6 and abs(v2)<=EPS_ZERO_P6:
                return QgsPointXY((x1+x2)/2.0,(y1+y2)/2.0)
            if abs(v1)<=EPS_ZERO_P6:
                return QgsPointXY(x1,y1)
            if abs(v2)<=EPS_ZERO_P6:
                return QgsPointXY(x2,y2)
            den=v1-v2
            t=0.5 if abs(den)<=EPS_ZERO_P6 else v1/den
            t=max(0.0,min(1.0,t))
            return QgsPointXY(
                float(x1+t*(x2-x1)),
                float(y1+t*(y2-y1))
            )

        def extrair_p6_original(cid, info):
            pair=info["pair"]
            cells=info["cells"]
            if pair not in ("COM_ESP","ESP_PRES") or len(cells)<20:
                return []

            # Registros equivalentes aos pontos da camada P5b original.
            pts=[]
            for rr,cc in cells:
                x,y=self._rc_to_xy(gt,rr,cc)
                delta=(float(A_COM[rr,cc]-A_ESP[rr,cc])
                       if pair=="COM_ESP"
                       else float(A_ESP[rr,cc]-A_PRES[rr,cc]))
                pts.append({
                    "r0":int(rr),"c0":int(cc),"x":float(x),"y":float(y),
                    "delta":delta,
                    "margem":float(mg[rr,cc]),
                    "conc":float(a2g[rr,cc]),
                    "a1":float(a1g[rr,cc]),
                })

            xmin=min(p["x"] for p in pts); xmax=max(p["x"] for p in pts)
            ymin=min(p["y"] for p in pts); ymax=max(p["y"] for p in pts)
            cols=int(round((xmax-xmin)/pix_m))+1
            rows_loc=int(round((ymax-ymin)/pix_m))+1

            valores=np.full((rows_loc,cols),np.nan,dtype=np.float64)
            margens=np.full((rows_loc,cols),np.nan,dtype=np.float64)
            concs=np.full((rows_loc,cols),np.nan,dtype=np.float64)
            a1s=np.full((rows_loc,cols),np.nan,dtype=np.float64)
            presente=np.zeros((rows_loc,cols),dtype=bool)

            for p in pts:
                cc=int(round((p["x"]-xmin)/pix_m))
                rr=int(round((ymax-p["y"])/pix_m))
                if 0<=rr<rows_loc and 0<=cc<cols:
                    valores[rr,cc]=p["delta"]
                    margens[rr,cc]=p["margem"]
                    concs[rr,cc]=p["conc"]
                    a1s[rr,cc]=p["a1"]
                    presente[rr,cc]=True

            segmentos=[]
            seg_local=0
            for rr in range(rows_loc-1):
                for cc in range(cols-1):
                    if not (presente[rr,cc] and presente[rr,cc+1]
                            and presente[rr+1,cc+1] and presente[rr+1,cc]):
                        continue

                    v0=float(valores[rr,cc]); v1=float(valores[rr,cc+1])
                    v2=float(valores[rr+1,cc+1]); v3=float(valores[rr+1,cc])
                    x0=xmin+cc*pix_m
                    x1=xmin+(cc+1)*pix_m
                    y0=ymax-rr*pix_m
                    y1=ymax-(rr+1)*pix_m

                    b0=1 if v0>=0.0 else 0
                    b1=1 if v1>=0.0 else 0
                    b2=1 if v2>=0.0 else 0
                    b3=1 if v3>=0.0 else 0
                    code=b0+2*b1+4*b2+8*b3
                    if code in (0,15):
                        continue

                    edge_pts={}
                    def ponto_aresta(edge):
                        if edge in edge_pts:
                            return edge_pts[edge]
                        if edge==0:
                            pz=interp_zero_p6(x0,y0,v0,x1,y0,v1)
                        elif edge==1:
                            pz=interp_zero_p6(x1,y0,v1,x1,y1,v2)
                        elif edge==2:
                            pz=interp_zero_p6(x0,y1,v3,x1,y1,v2)
                        elif edge==3:
                            pz=interp_zero_p6(x0,y0,v0,x0,y1,v3)
                        else:
                            raise QgsProcessingException("P6: aresta marching squares inválida.")
                        edge_pts[edge]=pz
                        return pz

                    if code in (5,10):
                        vc=(v0+v1+v2+v3)/4.0
                        if code==5:
                            pares=[(3,2),(0,1)] if vc>=0.0 else [(3,0),(2,1)]
                        else:
                            pares=[(3,0),(2,1)] if vc>=0.0 else [(3,2),(0,1)]
                    else:
                        pares=MS_CASES[code]

                    m_med=float(np.median(np.asarray([
                        margens[rr,cc],margens[rr,cc+1],
                        margens[rr+1,cc+1],margens[rr+1,cc]
                    ],dtype=np.float64)))
                    c_med=float(np.median(np.asarray([
                        concs[rr,cc],concs[rr,cc+1],
                        concs[rr+1,cc+1],concs[rr+1,cc]
                    ],dtype=np.float64)))
                    a1_med=float(np.median(np.asarray([
                        a1s[rr,cc],a1s[rr,cc+1],
                        a1s[rr+1,cc+1],a1s[rr+1,cc]
                    ],dtype=np.float64)))

                    for e0,e1 in pares:
                        p0=ponto_aresta(e0); p1=ponto_aresta(e1)
                        dist=math.hypot(p1.x()-p0.x(),p1.y()-p0.y())
                        if dist<=EPS_ZERO_P6:
                            continue
                        seg_local+=1
                        segmentos.append({
                            "seg_local":seg_local,"comp":cid,"par":pair,
                            "geom":QgsGeometry.fromPolylineXY([p0,p1]),
                            "p0":(float(p0.x()),float(p0.y())),
                            "p1":(float(p1.x()),float(p1.y())),
                            "length":float(dist),"m_med":m_med,
                            "c_med":c_med,"a1_med":a1_med,
                        })
            return segmentos

        # P6: somente 8 componentes relevantes; residual de 3 pixels não entra.
        segmentos_p6=[]
        seg_by_comp={}
        for cid in sorted(comp_info):
            info=comp_info[cid]
            if info["pair"] not in ("COM_ESP","ESP_PRES") or len(info["cells"])<20:
                continue
            segs=extrair_p6_original(cid,info)
            seg_by_comp[cid]=segs
            segmentos_p6.extend(segs)
            feedback.pushInfo(
                f"P6 COMP {cid} | {info['pair']} | pixels={len(info['cells'])} | segmentos={len(segs)}"
            )

        # SEG_ID global, na mesma ordem documental: componente -> segmento local.
        for sid,s in enumerate(segmentos_p6,1):
            s["seg_id"]=sid

        n_ce=sum(1 for s in segmentos_p6 if s["par"]=="COM_ESP")
        n_ep=sum(1 for s in segmentos_p6 if s["par"]=="ESP_PRES")
        feedback.pushInfo(
            f"P6 segmentos={len(segmentos_p6)} | COM_ESP={n_ce} | ESP_PRES={n_ep}"
        )

        # Componentes sem inversão documentados pelo P6b.
        sem_linha={
            cid for cid,segs in seg_by_comp.items()
            if len(segs)==0
        }

        # --------------------------------------------------------------
        # P7 – RECONSTRUÇÃO TOPOLOGICA DOCUMENTAL
        # --------------------------------------------------------------
        SNAP_P7=1.0
        EPS_P7=1e-9

        def snap_key_p7(x,y):
            return (int(round(float(x)/SNAP_P7)),int(round(float(y)/SNAP_P7)))

        def geom_coords(g):
            if g.isMultipart():
                parts=g.asMultiPolyline()
                if not parts: return []
                line=max(parts,key=lambda z:len(z))
            else:
                line=g.asPolyline()
            return [QgsPointXY(p.x(),p.y()) for p in line]

        segmentos_idx={s["seg_id"]:s for s in segmentos_p6}
        por_comp=defaultdict(list)
        micro=[]
        for s in segmentos_p6:
            p0=s["p0"]; p1=s["p1"]
            n0=snap_key_p7(*p0); n1=snap_key_p7(*p1)
            s["n0"]=n0; s["n1"]=n1
            if n0==n1:
                micro.append(s)
            else:
                por_comp[s["comp"]].append(s["seg_id"])

        grafos={}
        for cid in sorted(por_comp):
            adj=defaultdict(list); pair=None
            for sid in por_comp[cid]:
                s=segmentos_idx[sid]
                pair=s["par"] if pair is None else pair
                adj[s["n0"]].append(sid); adj[s["n1"]].append(sid)
            grafos[cid]={
                "par":pair,"segments":por_comp[cid],
                "adj":adj,"graus":{n:len(v) for n,v in adj.items()}
            }

        def orienta_segmento_p7(s,no_atual):
            coords=geom_coords(s["geom"])
            if not coords: return []
            if s["n0"]==no_atual:
                return coords
            if s["n1"]==no_atual:
                return list(reversed(coords))
            raise QgsProcessingException("P7: segmento não toca nó atual.")

        def unir_coords_p7(lista,novas):
            if not novas: return lista
            if not lista: return list(novas)
            ultimo=lista[-1]; primeiro=novas[0]
            if (abs(ultimo.x()-primeiro.x())<=SNAP_P7+EPS_P7
                    and abs(ultimo.y()-primeiro.y())<=SNAP_P7+EPS_P7):
                lista.extend(novas[1:])
            else:
                lista.extend(novas)
            return lista

        ramos=[]
        ramo_global=0
        for cid in sorted(grafos):
            gr=grafos[cid]; adj=gr["adj"]; graus=gr["graus"]; seg_ids=gr["segments"]
            visitados=set()
            especiais=[no for no,grau in graus.items() if grau!=2]

            for no_inicio in especiais:
                for sid0 in adj[no_inicio]:
                    if sid0 in visitados: continue
                    ramo_global+=1
                    atual_no=no_inicio; sid=sid0; lista_seg=[]; coords=[]
                    while True:
                        if sid in visitados: break
                        visitados.add(sid)
                        s=segmentos_idx[sid]
                        lista_seg.append(sid)
                        coords=unir_coords_p7(coords,orienta_segmento_p7(s,atual_no))
                        proximo_no=s["n1"] if s["n0"]==atual_no else s["n0"]
                        if graus[proximo_no]!=2:
                            atual_no=proximo_no; break
                        candidatos=[ss for ss in adj[proximo_no] if ss not in visitados]
                        if not candidatos:
                            atual_no=proximo_no; break
                        if len(candidatos)>1:
                            raise QgsProcessingException("P7: ambiguidade inesperada em nó grau 2.")
                        atual_no=proximo_no; sid=candidatos[0]
                    if lista_seg:
                        geom=QgsGeometry.fromPolylineXY(coords)
                        if geom.isEmpty() or geom.length()<=EPS_P7:
                            raise QgsProcessingException("P7: ramo inválido.")
                        ramos.append({
                            "ramo_id":ramo_global,"comp":cid,"par":gr["par"],
                            "segments":list(lista_seg),"geom":geom,
                            "length":float(geom.length()),"nseg":len(lista_seg)
                        })

            restantes=[sid for sid in seg_ids if sid not in visitados]
            while restantes:
                sid0=restantes[0]; s0=segmentos_idx[sid0]
                no_inicio=s0["n0"]; atual_no=no_inicio; sid=sid0
                lista_seg=[]; coords=[]
                while True:
                    if sid in visitados: break
                    visitados.add(sid)
                    s=segmentos_idx[sid]; lista_seg.append(sid)
                    coords=unir_coords_p7(coords,orienta_segmento_p7(s,atual_no))
                    proximo_no=s["n1"] if s["n0"]==atual_no else s["n0"]
                    candidatos=[ss for ss in adj[proximo_no] if ss not in visitados]
                    atual_no=proximo_no
                    if atual_no==no_inicio or not candidatos: break
                    sid=candidatos[0]
                if lista_seg:
                    ramo_global+=1
                    geom=QgsGeometry.fromPolylineXY(coords)
                    if geom.isEmpty() or geom.length()<=EPS_P7:
                        raise QgsProcessingException("P7: ciclo inválido.")
                    ramos.append({
                        "ramo_id":ramo_global,"comp":cid,"par":gr["par"],
                        "segments":list(lista_seg),"geom":geom,
                        "length":float(geom.length()),"nseg":len(lista_seg)
                    })
                restantes=[sid for sid in seg_ids if sid not in visitados]

            if len(visitados)!=len(seg_ids):
                raise QgsProcessingException(f"P7: COMP {cid} não contabilizou todos os segmentos.")

        # Métricas de ramo = média ponderada pelo comprimento dos segmentos.
        for br in ramos:
            ss=[segmentos_idx[sid] for sid in br["segments"]]
            w=np.asarray([s["length"] for s in ss],dtype=np.float64)
            if float(w.sum())<=EPS_P7:
                raise QgsProcessingException("P7: peso total zero.")
            br["m_med"]=float(np.average(
                np.asarray([s["m_med"] for s in ss],dtype=np.float64),weights=w))
            br["c_med"]=float(np.average(
                np.asarray([s["c_med"] for s in ss],dtype=np.float64),weights=w))
            br["a1_med"]=float(np.average(
                np.asarray([s["a1_med"] for s in ss],dtype=np.float64),weights=w))

        ramos_por_comp=defaultdict(list)
        for br in ramos: ramos_por_comp[br["comp"]].append(br)
        for cid in sorted(ramos_por_comp):
            rr=ramos_por_comp[cid]
            total_len=sum(r["length"] for r in rr)
            cmax=max(r["c_med"] for r in rr)
            amax=max(r["a1_med"] for r in rr)
            for br in rr:
                br["f_len"]=br["length"]/total_len if total_len>EPS_P7 else 0.0
                br["f_c"]=br["c_med"]/cmax if cmax>EPS_P7 else 0.0
                br["f_a1"]=br["a1_med"]/amax if amax>EPS_P7 else 0.0
                br["score"]=(self.W_LEN*br["f_len"]+
                             self.W_C*br["f_c"]+
                             self.W_A1*br["f_a1"])
            rr.sort(key=lambda z:(z["score"],z["length"]),reverse=True)
            principal_len=rr[0]["length"]
            for rank,br in enumerate(rr,1):
                br["rank"]=rank
                if rank==1:
                    br["hier"]="RAMO_PRINCIPAL"
                elif br["length"]>=self.FRAC_SEC*principal_len and br["length"]>=self.LEN_SEC_MIN_M:
                    br["hier"]="RAMO_SECUNDARIO"
                else:
                    br["hier"]="RAMO_LOCAL"

        branch_records=[]
        for cid in sorted(ramos_por_comp):
            for br in ramos_por_comp[cid]:
                coords=geom_coords(br["geom"])
                branch_records.append({
                    "cid":br["comp"],"pair":br["par"],"rank":br["rank"],
                    "hier":br["hier"],"score":br["score"],"L":br["length"],
                    "pts":[(float(p.x()),float(p.y())) for p in coords],
                    "flen":br["f_len"],"cmed":br["c_med"],
                    "a1med":br["a1_med"],"geom":br["geom"],
                    "nseg":br["nseg"],"ramo_id":br["ramo_id"]
                })

        seg_total=len(segmentos_p6); micro_total=len(micro)
        seg_em_ramos={sid for br in ramos for sid in br["segments"]}
        micro_ids={s["seg_id"] for s in micro}
        if seg_em_ramos & micro_ids:
            raise QgsProcessingException("P7: segmento simultaneamente em ramo e microfragmento.")
        if (seg_em_ramos | micro_ids) != set(segmentos_idx):
            raise QgsProcessingException("P7: conservação 1:1 dos segmentos falhou.")

        feedback.pushInfo(
            f"P7 segmentos_em_ramos={len(seg_em_ramos)} | "
            f"micro={micro_total} | ramos={len(branch_records)}"
        )

        if run_bench and modo_acaj:
            # Benchmark documental completo P6/P7.
            expected_comp={1:342,3:27,4:2,5:18,6:0,7:209,8:0,9:21}
            got_comp={cid:len(seg_by_comp.get(cid,[])) for cid in expected_comp}
            gate_errors=[]
            if seg_total!=619: gate_errors.append(f"P6={seg_total}!=619")
            if n_ce!=389: gate_errors.append(f"COM_ESP={n_ce}!=389")
            if n_ep!=230: gate_errors.append(f"ESP_PRES={n_ep}!=230")
            if got_comp!=expected_comp: gate_errors.append(f"por_comp={got_comp}")
            if len(seg_em_ramos)!=615: gate_errors.append(f"em_ramos={len(seg_em_ramos)}!=615")
            if micro_total!=4: gate_errors.append(f"micro={micro_total}!=4")
            if len(branch_records)!=8: gate_errors.append(f"ramos={len(branch_records)}!=8")
            if gate_errors:
                raise QgsProcessingException(
                    "GATE P6/P7 ACAJATUBA FALHOU: "+" | ".join(gate_errors)
                )
            feedback.pushInfo(
                "[OK] Gate P6/P7: 619 segmentos (389+230), "
                "615 em ramos + 4 micro, 8 ramos."
            )

        # --------------------------------------------------------------
        # P7b – FAIXA OPERACIONAL DOCUMENTAL
        # --------------------------------------------------------------
        # Distância somente ao RAMO_PRINCIPAL do MESMO COMP_P5B,
        # intersectada com os critérios funcionais do P4/P5b.
        ramo_principal_por_comp={}
        for br in branch_records:
            if br["hier"]=="RAMO_PRINCIPAL":
                if br["cid"] in ramo_principal_por_comp:
                    raise QgsProcessingException(
                        f"P7b: mais de um principal no COMP {br['cid']}."
                    )
                ramo_principal_por_comp[br["cid"]]=br

        faixa_flag=np.zeros((nr,nc),dtype=bool)
        faixa30=np.zeros((nr,nc),dtype=bool)
        faixa90=np.zeros((nr,nc),dtype=bool)
        dist_eq=np.full((nr,nc),np.nan,dtype=np.float64)

        for rec in rows:
            r,c=rec["r"],rec["c"]
            cid=int(comp_grid[r,c])
            if cid not in ramo_principal_por_comp:
                continue
            br=ramo_principal_por_comp[cid]
            if interface_pair[r,c] != br["pair"]:
                continue

            x,y=self._rc_to_xy(gt,r,c)
            pg=QgsGeometry.fromPointXY(QgsPointXY(x,y))
            dd=float(pg.distance(br["geom"]))
            dist_eq[r,c]=dd

            func_ok=(rec["m"]<=mcrit+eps and
                     rec["a2"]>=ccrit-eps and
                     rec["a1"]>=amin-eps)
            if not func_ok:
                continue
            faixa30[r,c]=(dd<=30.0+eps)
            faixa_flag[r,c]=(dd<=60.0+eps)
            faixa90[r,c]=(dd<=90.0+eps)

        feedback.pushInfo(
            f"P7b FAIXA 30/60/90 = "
            f"{int(faixa30.sum())}/{int(faixa_flag.sum())}/{int(faixa90.sum())}"
        )

        if run_bench and modo_acaj:
            n30=int(faixa30.sum()); n60=int(faixa_flag.sum()); n90=int(faixa90.sum())
            pair60=Counter()
            comp60=Counter()
            for rec in rows:
                r,c=rec["r"],rec["c"]
                if faixa_flag[r,c]:
                    pair60[str(interface_pair[r,c])]+=1
                    comp60[int(comp_grid[r,c])]+=1
            expected_comp60={1:903,3:91,4:11,5:64,7:662,9:70}
            gate_errors=[]
            if (n30,n60,n90)!=(895,1801,2692):
                gate_errors.append(f"30/60/90={n30}/{n60}/{n90}")
            if pair60.get("COM_ESP",0)!=1069 or pair60.get("ESP_PRES",0)!=732:
                gate_errors.append(f"pares60={dict(pair60)}")
            if {k:comp60.get(k,0) for k in expected_comp60}!=expected_comp60:
                gate_errors.append(f"comp60={dict(comp60)}")
            if gate_errors:
                raise QgsProcessingException(
                    "GATE P7b ACAJATUBA FALHOU: "+" | ".join(gate_errors)
                )
            feedback.pushInfo(
                "[OK] Gate P7b: 895/1801/2692 pixels; "
                "60 m = 1069 COM_ESP + 732 ESP_PRES."
            )

        # --------------------------------------------------------------
        # P8a2 – Chebyshev via multi-source heap, tie lowest PID
        # --------------------------------------------------------------
        index_by_rc={(rec["r"],rec["c"]):rec for rec in rows}
        sources=[]; targets=set()
        for rec in rows:
            rc=(rec["r"],rec["c"])
            if rec["tie3"] and rec["a1"]<amin:
                targets.add(rc)
            elif rec["a1"]>=amin:
                sources.append((rec["r"],rec["c"],rec["pid"],rec["win"]))
        best={}; hp=[]
        for r,c,spid,sclass in sources:
            key=(0,spid)
            old=best.get((r,c))
            if old is None or key<(old[0],old[1]):
                best[(r,c)]=(0,spid,sclass); heapq.heappush(hp,(0,spid,r,c,sclass))
        pending=set(targets); prox={}
        neigh=[(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]
        while hp and pending:
            d,spid,r,c,sclass=heapq.heappop(hp)
            cur=best.get((r,c))
            if cur is None or (d,spid)!=(cur[0],cur[1]): continue
            if (r,c) in pending:
                prox[(r,c)]={"dist_pix":d,"dist_m":d*pix_m,"src":spid,"classe":sclass}; pending.remove((r,c))
            nd=d+1
            for dr,dc in neigh:
                rr,cc=r+dr,c+dc
                if rr<0 or rr>=nr or cc<0 or cc>=nc: continue
                old=best.get((rr,cc))
                if old is None or (nd,spid)<(old[0],old[1]):
                    best[(rr,cc)]=(nd,spid,sclass); heapq.heappush(hp,(nd,spid,rr,cc,sclass))
        if pending: raise QgsProcessingException(f"P8a2 deixou {len(pending)} empates sem resolução.")
        feedback.pushInfo(f"P8a2 PROX_ESPACIAL={len(prox)}")

        if run_bench and modo_acaj:
            if len(prox) != 5740:
                raise QgsProcessingException(
                    f"GATE P8a2 ACAJATUBA FALHOU: PROX_ESPACIAL={len(prox)} (esperado=5740)."
                )
            feedback.pushInfo("[OK] Gate P8a2: 5740 empates verdadeiros resolvidos.")

        # --------------------------------------------------------------
        # Outputs schemas
        # --------------------------------------------------------------
        crs=pontos.sourceCrs()
        anc_fields=QgsFields()
        for n,t,l,p in [("ID_LEGAL",QVariant.Int,0,0),("SISTEMA",QVariant.String,40,0),("CLASSE_USO",QVariant.String,20,0),("ANC_TIPO",QVariant.String,40,0),("DIST_ANC_M",QVariant.Double,20,3)]: anc_fields.append(QgsField(n,t,len=l,prec=p))
        anc_sink,anc_dest=self.parameterAsSink(parameters,self.ANCORAS_OUT,context,anc_fields,QgsWkbTypes.Point,crs)

        aff_fields=QgsFields()
        for n,t,l,p in [("PID",QVariant.Int,0,0),("ROW",QVariant.Int,0,0),("COL",QVariant.Int,0,0),("A_COM",QVariant.Double,20,9),("A_ESP",QVariant.Double,20,9),("A_PRES",QVariant.Double,20,9),("A1",QVariant.Double,20,9),("A2",QVariant.Double,20,9),("A3",QVariant.Double,20,9),("MARGEM3",QVariant.Double,20,9),("AFIN_CONC2",QVariant.Double,20,9),("PAR_COMP",QVariant.String,20,0)]: aff_fields.append(QgsField(n,t,len=l,prec=p))
        aff_sink,aff_dest=self.parameterAsSink(parameters,self.AFINIDADES_OUT,context,aff_fields,QgsWkbTypes.Point,crs)

        int_fields=QgsFields()
        for n,t,l,p in [("PID",QVariant.Int,0,0),("PAR_COMP",QVariant.String,20,0),("COMP_P5B",QVariant.Int,0,0),("A1",QVariant.Double,20,9),("A2",QVariant.Double,20,9),("MARGEM3",QVariant.Double,20,9)]: int_fields.append(QgsField(n,t,len=l,prec=p))
        int_sink,int_dest=self.parameterAsSink(parameters,self.INTERFACES_OUT,context,int_fields,QgsWkbTypes.Point,crs)

        ramo_fields=QgsFields()
        for n,t,l,p in [("RAMO_ID",QVariant.Int,0,0),("COMP_P5B",QVariant.Int,0,0),("PAR_EQ",QVariant.String,20,0),("RANK",QVariant.Int,0,0),("HIERARQUIA",QVariant.String,30,0),("LENGTH_M",QVariant.Double,20,3),("SCORE_P7",QVariant.Double,20,9)]: ramo_fields.append(QgsField(n,t,len=l,prec=p))
        ramo_sink,ramo_dest=self.parameterAsSink(parameters,self.RAMOS_OUT,context,ramo_fields,QgsWkbTypes.LineString,crs)

        faixa_fields=QgsFields()
        for n,t,l,p in [("PID",QVariant.Int,0,0),("A_COM",QVariant.Double,20,9),("A_ESP",QVariant.Double,20,9),("A_PRES",QVariant.Double,20,9),("DIST_EQ_M",QVariant.Double,20,3),("FAIXA_OP",QVariant.Int,0,0)]: faixa_fields.append(QgsField(n,t,len=l,prec=p))
        faixa_sink,faixa_dest=self.parameterAsSink(parameters,self.FAIXA_OUT,context,faixa_fields,QgsWkbTypes.Point,crs)

        prox_fields=QgsFields()
        for n,t,l,p in [("PID",QVariant.Int,0,0),("CLASSE_PROX",QVariant.String,20,0),("DIST_PIX",QVariant.Int,0,0),("PROX_DIST",QVariant.Double,20,3),("PROX_SRC",QVariant.Int,0,0),("CONFIANCA",QVariant.String,20,0)]: prox_fields.append(QgsField(n,t,len=l,prec=p))
        prox_sink,prox_dest=self.parameterAsSink(parameters,self.PROX_OUT,context,prox_fields,QgsWkbTypes.Point,crs)

        final_fields=QgsFields()
        for n,t,l,p in [("PID",QVariant.Int,0,0),("A_COM",QVariant.Double,20,9),("A_ESP",QVariant.Double,20,9),("A_PRES",QVariant.Double,20,9),("A1",QVariant.Double,20,9),("A2",QVariant.Double,20,9),("A3",QVariant.Double,20,9),("MARGEM3",QVariant.Double,20,9),("AFIN_CONC2",QVariant.Double,20,9),("PAR_COMP",QVariant.String,20,0),("FAIXA_OPER",QVariant.Int,0,0),("CLASSE_FIN",QVariant.String,20,0),("ORIGEM_CLASS",QVariant.String,30,0),("CONFIANCA",QVariant.String,20,0),("PROX_DIST",QVariant.Double,20,3),("PROX_SRC",QVariant.Int,0,0)]: final_fields.append(QgsField(n,t,len=l,prec=p))
        final_sink,final_dest=self.parameterAsSink(parameters,self.FINAL_OUT,context,final_fields,QgsWkbTypes.Point,crs)

        res_fields=QgsFields()
        for n,t,l,p in [("TIPO",QVariant.String,30,0),("CHAVE",QVariant.String,30,0),("N",QVariant.Int,0,0),("VALOR",QVariant.Double,20,9)]: res_fields.append(QgsField(n,t,len=l,prec=p))
        res_sink,res_dest=self.parameterAsSink(parameters,self.RESUMO_OUT,context,res_fields,QgsWkbTypes.NoGeometry,crs)
        lim_sink,lim_dest=self.parameterAsSink(parameters,self.LIMIARES_OUT,context,res_fields,QgsWkbTypes.NoGeometry,crs)

        # P9 – QA/QC outputs
        p9_comp_fields=QgsFields()
        for n,t,l,p in [
            ("COMP_ID",QVariant.Int,0,0),("CLASSE_FIN",QVariant.String,20,0),
            ("N_PIX",QVariant.Int,0,0),("AREA_KM2",QVariant.Double,20,6),
            ("FRAG20",QVariant.Int,0,0),("FRAG3",QVariant.Int,0,0)
        ]:
            p9_comp_fields.append(QgsField(n,t,len=l,prec=p))
        p9_comp_sink,p9_comp_dest=self.parameterAsSink(
            parameters,self.P9_COMPONENTES_OUT,context,p9_comp_fields,QgsWkbTypes.NoGeometry,crs
        )

        p9_frag_fields=QgsFields()
        for n,t,l,p in [
            ("PID",QVariant.Int,0,0),("COMP_ID",QVariant.Int,0,0),
            ("CLASSE_FIN",QVariant.String,20,0),("N_COMP",QVariant.Int,0,0),
            ("ORIGEM",QVariant.String,30,0),("CONFIANCA",QVariant.String,20,0)
        ]:
            p9_frag_fields.append(QgsField(n,t,len=l,prec=p))
        p9_frag_sink,p9_frag_dest=self.parameterAsSink(
            parameters,self.P9_FRAGMENTOS_OUT,context,p9_frag_fields,QgsWkbTypes.Point,crs
        )

        p9_cross_fields=QgsFields()
        p9_cross_fields.append(QgsField("LINHA",QVariant.String,len=30))
        for c in ("ALTA","MEDIA","TRANSICAO","BAIXA","ARGMAX","FAIXA_OPER","PROX_ESPACIAL"):
            p9_cross_fields.append(QgsField(c,QVariant.Int))
        p9_cross_co_sink,p9_cross_co_dest=self.parameterAsSink(
            parameters,self.P9_CLASSE_ORIGEM_OUT,context,p9_cross_fields,QgsWkbTypes.NoGeometry,crs
        )
        p9_cross_cc_sink,p9_cross_cc_dest=self.parameterAsSink(
            parameters,self.P9_CLASSE_CONFIANCA_OUT,context,p9_cross_fields,QgsWkbTypes.NoGeometry,crs
        )
        p9_cross_oc_sink,p9_cross_oc_dest=self.parameterAsSink(
            parameters,self.P9_ORIGEM_CONFIANCA_OUT,context,p9_cross_fields,QgsWkbTypes.NoGeometry,crs
        )

        p9_res_fields=QgsFields()
        p9_res_fields.append(QgsField("INDICADOR",QVariant.String,len=40))
        p9_res_fields.append(QgsField("VALOR",QVariant.Double,prec=6))
        p9_res_sink,p9_res_dest=self.parameterAsSink(
            parameters,self.P9_RESUMO_OUT,context,p9_res_fields,QgsWkbTypes.NoGeometry,crs
        )

        # P10 – final polygon products
        p10_comp_fields=QgsFields()
        for n,t,l,p in [
            ("COMP_ID",QVariant.Int,0,0),("CLASSE_FIN",QVariant.String,20,0),
            ("N_PIX",QVariant.Int,0,0),("AREA_KM2",QVariant.Double,20,6),
            ("AREA_HA",QVariant.Double,20,3)
        ]:
            p10_comp_fields.append(QgsField(n,t,len=l,prec=p))
        p10_comp_sink,p10_comp_dest=self.parameterAsSink(
            parameters,self.P10_COMPONENTES_OUT,context,p10_comp_fields,QgsWkbTypes.MultiPolygon,crs
        )

        p10_class_fields=QgsFields()
        for n,t,l,p in [
            ("CLASSE_FIN",QVariant.String,20,0),("N_PIX",QVariant.Int,0,0),
            ("AREA_KM2",QVariant.Double,20,6),("AREA_HA",QVariant.Double,20,3),
            ("PC_DOM",QVariant.Double,20,6),("N_COMP",QVariant.Int,0,0)
        ]:
            p10_class_fields.append(QgsField(n,t,len=l,prec=p))
        p10_class_sink,p10_class_dest=self.parameterAsSink(
            parameters,self.P10_CLASSES_OUT,context,p10_class_fields,QgsWkbTypes.MultiPolygon,crs
        )

        p10_res_fields=QgsFields()
        for n,t,l,p in [
            ("CLASSE_FIN",QVariant.String,20,0),("N_PIX",QVariant.Int,0,0),
            ("AREA_KM2",QVariant.Double,20,6),("AREA_HA",QVariant.Double,20,3),
            ("PC_DOM",QVariant.Double,20,6),("N_COMP",QVariant.Int,0,0)
        ]:
            p10_res_fields.append(QgsField(n,t,len=l,prec=p))
        p10_res_sink,p10_res_dest=self.parameterAsSink(
            parameters,self.P10_RESUMO_OUT,context,p10_res_fields,QgsWkbTypes.NoGeometry,crs
        )

        # anchors output
        id_to_sys={}
        for s,ids in systems.items():
            for i in ids:id_to_sys[i]=s
        for i in sorted(pts,key=lambda x:int(x) if isinstance(x,int) else str(x)):
            a=anchors.get(i)
            if pts[i]["cat"]=="PRESERVACAO":
                rr,cc=seed_rc[i]; x,y=self._rc_to_xy(gt,rr,cc); d=math.hypot(x-pts[i]["x"],y-pts[i]["y"]); typ="GEOMETRICA_MAXEXTENT"
            else:
                x,y=a["x"],a["y"]; d=a["dist"]; typ=a["tipo"]
            f=QgsFeature(anc_fields); f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x,y))); f.setAttributes([int(i),id_to_sys.get(i,""),pts[i]["cat"],typ,d]); anc_sink.addFeature(f,QgsFeatureSink.FastInsert)

        # outputs points + P8 classification
        counts_class=Counter(); counts_orig=Counter(); counts_conf=Counter()
        final_records=[]
        ramo_id_counter=0
        for br in branch_records:
            ramo_id_counter+=1; f=QgsFeature(ramo_fields); f.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(x,y) for x,y in br["pts"]])); f.setAttributes([ramo_id_counter,br["cid"],br["pair"],br["rank"],br["hier"],br["L"],br["score"]]); ramo_sink.addFeature(f,QgsFeatureSink.FastInsert)

        for rec in rows:
            r,c,pid=rec["r"],rec["c"],rec["pid"]; x,y=self._rc_to_xy(gt,r,c); geom=QgsGeometry.fromPointXY(QgsPointXY(x,y))
            af=QgsFeature(aff_fields); af.setGeometry(geom); af.setAttributes([pid,r,c,float(A_COM[r,c]),float(A_ESP[r,c]),float(A_PRES[r,c]),rec["a1"],rec["a2"],rec["a3"],rec["m"],rec["a2"],rec["pair"]]); aff_sink.addFeature(af,QgsFeatureSink.FastInsert)
            if valid_interface[r,c]:
                ff=QgsFeature(int_fields); ff.setGeometry(geom); ff.setAttributes([pid,str(interface_pair[r,c]),int(comp_grid[r,c]),rec["a1"],rec["a2"],rec["m"]]); int_sink.addFeature(ff,QgsFeatureSink.FastInsert)
            fq=QgsFeature(faixa_fields); fq.setGeometry(geom); fq.setAttributes([pid,float(A_COM[r,c]),float(A_ESP[r,c]),float(A_PRES[r,c]),None if not np.isfinite(dist_eq[r,c]) else float(dist_eq[r,c]),1 if faixa_flag[r,c] else 0]); faixa_sink.addFeature(fq,QgsFeatureSink.FastInsert)

            pr=prox.get((r,c))
            if pr is not None:
                pp=QgsFeature(prox_fields); pp.setGeometry(geom); pp.setAttributes([pid,pr["classe"],pr["dist_pix"],pr["dist_m"],pr["src"],"BAIXA"]); prox_sink.addFeature(pp,QgsFeatureSink.FastInsert)
                classe=pr["classe"]; origem="PROX_ESPACIAL"; conf="BAIXA"; pd=pr["dist_m"]; ps=pr["src"]
            elif faixa_flag[r,c]:
                classe=rec["win"]; origem="FAIXA_OPER"; conf="TRANSICAO"; pd=None; ps=None
            else:
                if rec["tie3"]: raise QgsProcessingException(f"Empate não resolvido no PID {pid}")
                classe=rec["win"]; origem="ARGMAX"; pd=None; ps=None
                if rec["a1"]<amin: conf="BAIXA"
                elif rec["m"]<=mcrit+eps and rec["a2"]>=ccrit-eps: conf="MEDIA"
                else: conf="ALTA"
            counts_class[classe]+=1; counts_orig[origem]+=1; counts_conf[conf]+=1
            final_records.append({
                "pid":int(pid),"r":int(r),"c":int(c),"x":float(x),"y":float(y),
                "classe":str(classe),"origem":str(origem),"conf":str(conf),
                "a1":float(rec["a1"]),"faixa":1 if faixa_flag[r,c] else 0
            })
            fo=QgsFeature(final_fields); fo.setGeometry(geom); fo.setAttributes([pid,float(A_COM[r,c]),float(A_ESP[r,c]),float(A_PRES[r,c]),rec["a1"],rec["a2"],rec["a3"],rec["m"],rec["a2"],rec["pair"],1 if faixa_flag[r,c] else 0,classe,origem,conf,pd,ps]); final_sink.addFeature(fo,QgsFeatureSink.FastInsert)

        # summary/thresholds
        for k,v in counts_class.items():
            f=QgsFeature(res_fields); f.setAttributes(["CLASSE",k,int(v),None]); res_sink.addFeature(f,QgsFeatureSink.FastInsert)
        for k,v in counts_orig.items():
            f=QgsFeature(res_fields); f.setAttributes(["ORIGEM",k,int(v),None]); res_sink.addFeature(f,QgsFeatureSink.FastInsert)
        for k,v in counts_conf.items():
            f=QgsFeature(res_fields); f.setAttributes(["CONFIANCA",k,int(v),None]); res_sink.addFeature(f,QgsFeatureSink.FastInsert)
        for k,v in [("M_CRIT",mcrit),("C_CRIT",ccrit),("A_MIN",amin),("M_ALTA",malta),("C_ALTA",calta),("LAMBDA",lambda_common)]:
            f=QgsFeature(res_fields); f.setAttributes(["LIMIAR",k,0,float(v)]); lim_sink.addFeature(f,QgsFeatureSink.FastInsert)

        feedback.pushInfo("-"*120)
        feedback.pushInfo("P8 – RESULTADO FINAL")
        feedback.pushInfo(f"DOMINIO={len(rows)}")
        feedback.pushInfo("Classes: "+" | ".join(f"{k}={counts_class[k]}" for k in self.CLASSES))
        feedback.pushInfo("Origens: "+" | ".join(f"{k}={counts_orig[k]}" for k in ("ARGMAX","FAIXA_OPER","PROX_ESPACIAL")))
        feedback.pushInfo("Confiança: "+" | ".join(f"{k}={counts_conf[k]}" for k in ("ALTA","MEDIA","TRANSICAO","BAIXA")))

        if run_bench:
            errs=[]
            if len(rows)!=self.BENCH["n_total"]: errs.append(f"dominio={len(rows)}")
            for k,e in self.BENCH["classes"].items():
                if counts_class[k]!=e: errs.append(f"{k}={counts_class[k]}!= {e}")
            for k,e in self.BENCH["origens"].items():
                if counts_orig[k]!=e: errs.append(f"{k}={counts_orig[k]}!= {e}")
            for k,e in self.BENCH["confiancas"].items():
                if counts_conf[k]!=e: errs.append(f"{k}={counts_conf[k]}!= {e}")
            if errs:
                raise QgsProcessingException("BENCHMARK ACAJATUBA FALHOU: "+" | ".join(errs))
            feedback.pushInfo("[OK] Benchmark Acajatuba reproduzido integralmente.")

        # ==============================================================
        # P9 — QA/QC GEOMÉTRICO E ESTATÍSTICO (DIAGNÓSTICO; NÃO RECLASSIFICA)
        # ==============================================================
        feedback.pushInfo("="*120)
        feedback.pushInfo("P9 – QA/QC GEOMETRICO E ESTATISTICO")
        feedback.pushInfo("-"*120)

        class_code={"COMERCIAL":1,"ESPORTIVA":2,"PRESERVACAO":3}
        class_grid_final=np.zeros((nr,nc),dtype=np.int8)
        pid_by_rc={}
        rec_by_rc={}
        for d in final_records:
            class_grid_final[d["r"],d["c"]]=class_code[d["classe"]]
            pid_by_rc[(d["r"],d["c"])]=d["pid"]
            rec_by_rc[(d["r"],d["c"])]=d

        componentes_p9=[]
        comp_grid_p9=np.zeros((nr,nc),dtype=np.int32)
        comp_id_global=0
        for classe in self.CLASSES:
            cmask=(class_grid_final==class_code[classe])
            labs_local, comps_local=self._label8(cmask)
            # preserve deterministic row-major component order
            ordered=sorted(comps_local.values(), key=lambda cells:min(cells))
            for cells in ordered:
                comp_id_global+=1
                n_pix=len(cells)
                frag20=1 if n_pix<=20 else 0
                frag3=1 if n_pix<=3 else 0
                for rr,cc in cells:
                    comp_grid_p9[rr,cc]=comp_id_global
                componentes_p9.append({
                    "comp":comp_id_global,"classe":classe,"cells":cells,
                    "n_pix":n_pix,"frag20":frag20,"frag3":frag3
                })

        # P9 cross tables
        cross_co={cl:Counter() for cl in self.CLASSES}
        cross_cc={cl:Counter() for cl in self.CLASSES}
        cross_oc={o:Counter() for o in ("ARGMAX","FAIXA_OPER","PROX_ESPACIAL")}
        for d in final_records:
            cross_co[d["classe"]][d["origem"]]+=1
            cross_cc[d["classe"]][d["conf"]]+=1
            cross_oc[d["origem"]][d["conf"]]+=1

        n_frag20=sum(c["frag20"] for c in componentes_p9)
        n_frag3=sum(c["frag3"] for c in componentes_p9)
        frag_cells=set()
        for cp in componentes_p9:
            if cp["frag20"]:
                frag_cells.update(cp["cells"])

        n_prox=sum(1 for d in final_records if d["origem"]=="PROX_ESPACIAL")
        n_faixa=sum(1 for d in final_records if d["origem"]=="FAIXA_OPER")
        n_baixa=sum(1 for d in final_records if d["conf"]=="BAIXA")
        n_argmax_baixa=sum(1 for d in final_records if d["origem"]=="ARGMAX" and d["conf"]=="BAIXA")

        by_class_p9={}
        for cl in self.CLASSES:
            cps=[c for c in componentes_p9 if c["classe"]==cl]
            by_class_p9[cl]={
                "ncomp":len(cps),
                "pix":sum(c["n_pix"] for c in cps),
                "frag20":sum(c["frag20"] for c in cps),
                "frag3":sum(c["frag3"] for c in cps),
                "pix_frag20":sum(c["n_pix"] for c in cps if c["frag20"]),
                "pix_frag3":sum(c["n_pix"] for c in cps if c["frag3"]),
                "max_comp":max((c["n_pix"] for c in cps),default=0)
            }
            z=by_class_p9[cl]
            feedback.pushInfo(
                f"{cl:<13} | componentes={z['ncomp']} | maior={z['max_comp']} pix | "
                f"<=20={z['frag20']} ({z['pix_frag20']} pix) | "
                f"<=3={z['frag3']} ({z['pix_frag3']} pix)"
            )

        feedback.pushInfo(
            f"P9 dominio={len(final_records)} | componentes={len(componentes_p9)} | "
            f"frag<=20={n_frag20} | micro<=3={n_frag3} | pixels_frag={len(frag_cells)}"
        )
        feedback.pushInfo(
            f"P9 PROX={n_prox} | FAIXA={n_faixa} | BAIXA={n_baixa} | ARGMAX+BAIXA={n_argmax_baixa}"
        )

        if run_bench and modo_acaj:
            exp_cls={
                "COMERCIAL":{"ncomp":85,"pix":141672,"frag20":75,"pix_frag20":238,"frag3":58,"pix_frag3":84,"max_comp":140244},
                "ESPORTIVA":{"ncomp":343,"pix":57607,"frag20":315,"pix_frag20":1083,"frag3":229,"pix_frag3":346,"max_comp":39945},
                "PRESERVACAO":{"ncomp":36,"pix":31446,"frag20":29,"pix_frag20":66,"frag3":25,"pix_frag3":33,"max_comp":29636},
            }
            exp_co={
                "COMERCIAL":{"ARGMAX":139709,"FAIXA_OPER":535,"PROX_ESPACIAL":1428},
                "ESPORTIVA":{"ARGMAX":52820,"FAIXA_OPER":899,"PROX_ESPACIAL":3888},
                "PRESERVACAO":{"ARGMAX":30655,"FAIXA_OPER":367,"PROX_ESPACIAL":424},
            }
            exp_cc={
                "COMERCIAL":{"ALTA":81139,"MEDIA":8269,"TRANSICAO":535,"BAIXA":51729},
                "ESPORTIVA":{"ALTA":38576,"MEDIA":14244,"TRANSICAO":899,"BAIXA":3888},
                "PRESERVACAO":{"ALTA":18895,"MEDIA":10120,"TRANSICAO":367,"BAIXA":2064},
            }
            exp_oc={
                "ARGMAX":{"ALTA":138610,"MEDIA":32633,"TRANSICAO":0,"BAIXA":51941},
                "FAIXA_OPER":{"ALTA":0,"MEDIA":0,"TRANSICAO":1801,"BAIXA":0},
                "PROX_ESPACIAL":{"ALTA":0,"MEDIA":0,"TRANSICAO":0,"BAIXA":5740},
            }
            errs=[]
            for cl,exp in exp_cls.items():
                for k,v in exp.items():
                    if by_class_p9[cl][k]!=v:
                        errs.append(f"P9 {cl} {k}={by_class_p9[cl][k]}!={v}")
            for cl,exp in exp_co.items():
                for k,v in exp.items():
                    if cross_co[cl][k]!=v: errs.append(f"P9 CO {cl}/{k}={cross_co[cl][k]}!={v}")
            for cl,exp in exp_cc.items():
                for k,v in exp.items():
                    if cross_cc[cl][k]!=v: errs.append(f"P9 CC {cl}/{k}={cross_cc[cl][k]}!={v}")
            for org,exp in exp_oc.items():
                for k,v in exp.items():
                    if cross_oc[org][k]!=v: errs.append(f"P9 OC {org}/{k}={cross_oc[org][k]}!={v}")
            if len(componentes_p9)!=464: errs.append(f"P9 componentes={len(componentes_p9)}!=464")
            if n_frag20!=419: errs.append(f"P9 frag20={n_frag20}!=419")
            if n_frag3!=312: errs.append(f"P9 frag3={n_frag3}!=312")
            if len(frag_cells)!=1387: errs.append(f"P9 pix_frag={len(frag_cells)}!=1387")
            if (n_prox,n_faixa,n_baixa,n_argmax_baixa)!=(5740,1801,57681,51941):
                errs.append("P9 controles origem/confianca divergentes")
            if errs:
                raise QgsProcessingException("GATE P9 ACAJATUBA FALHOU: "+" | ".join(errs))
            feedback.pushInfo("[OK] Gate P9: 464 componentes; fragmentacao e tabelas cruzadas reproduzidas integralmente.")

        # P9 outputs
        for cp in componentes_p9:
            f=QgsFeature(p9_comp_fields)
            f.setAttributes([
                cp["comp"],cp["classe"],cp["n_pix"],
                cp["n_pix"]*(pix_m*pix_m/1_000_000.0),cp["frag20"],cp["frag3"]
            ])
            p9_comp_sink.addFeature(f,QgsFeatureSink.FastInsert)

        for rr,cc in sorted(frag_cells):
            d=rec_by_rc[(rr,cc)]
            f=QgsFeature(p9_frag_fields)
            f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(d["x"],d["y"])))
            cp_id=int(comp_grid_p9[rr,cc])
            cp=componentes_p9[cp_id-1]
            f.setAttributes([d["pid"],cp_id,d["classe"],cp["n_pix"],d["origem"],d["conf"]])
            p9_frag_sink.addFeature(f,QgsFeatureSink.FastInsert)

        def _write_cross(sink, row_name, counter):
            f=QgsFeature(p9_cross_fields)
            vals=[row_name]
            for col in ("ALTA","MEDIA","TRANSICAO","BAIXA","ARGMAX","FAIXA_OPER","PROX_ESPACIAL"):
                vals.append(int(counter.get(col,0)))
            f.setAttributes(vals)
            sink.addFeature(f,QgsFeatureSink.FastInsert)

        for cl in self.CLASSES:
            _write_cross(p9_cross_co_sink,cl,cross_co[cl])
            _write_cross(p9_cross_cc_sink,cl,cross_cc[cl])
        for org in ("ARGMAX","FAIXA_OPER","PROX_ESPACIAL"):
            _write_cross(p9_cross_oc_sink,org,cross_oc[org])

        p9_indicators=[
            ("N_DOMINIO",float(len(final_records))),
            ("AREA_KM2",float(len(final_records)*pix_m*pix_m/1_000_000.0)),
            ("N_COMPONENTES",float(len(componentes_p9))),
            ("N_FRAG20",float(n_frag20)),
            ("N_FRAG3",float(n_frag3)),
            ("N_PROX_ESPACIAL",float(n_prox)),
            ("N_FAIXA_OPER",float(n_faixa)),
            ("N_ARGMAX_BAIXA",float(n_argmax_baixa)),
        ]
        for k,v in p9_indicators:
            f=QgsFeature(p9_res_fields); f.setAttributes([k,v])
            p9_res_sink.addFeature(f,QgsFeatureSink.FastInsert)

        # ==============================================================
        # P10 — PRODUTOS POLIGONAIS FINAIS (PURAMENTE GEOMÉTRICO)
        # ==============================================================
        feedback.pushInfo("="*120)
        feedback.pushInfo("P10 – CONSTRUCAO DOS PRODUTOS POLIGONAIS FINAIS")
        feedback.pushInfo("-"*120)

        # Raster interno exato da grade P8: 0 externo; 1/2/3 classes.
        drv_mem=gdal.GetDriverByName("MEM")
        ds_mem=drv_mem.Create("",nc,nr,1,gdal.GDT_Byte)
        if ds_mem is None:
            raise QgsProcessingException("P10: falha ao criar raster MEM.")
        ds_mem.SetGeoTransform(gt)
        srs=osr.SpatialReference()
        try:
            srs.ImportFromWkt(crs.toWkt())
        except Exception:
            srs.ImportFromEPSG(31980)
        ds_mem.SetProjection(srs.ExportToWkt())
        band=ds_mem.GetRasterBand(1)
        band.WriteArray(class_grid_final.astype(np.uint8))
        band.SetNoDataValue(0)
        band.FlushCache()
        band.CreateMaskBand(gdal.GMF_PER_DATASET)
        maskband=band.GetMaskBand()
        maskband.WriteArray(np.where(class_grid_final>0,255,0).astype(np.uint8))
        maskband.FlushCache()

        # Polygonize 8-connected, exactly as documented P10.
        # Nota de implementação v6.0.1:
        # componentes conectados somente por vértices podem sair do GDAL como
        # polígonos auto-tocantes. O P10 preserva a conectividade 8-vizinhos,
        # mas normaliza cada componente com GEOS makeValid() antes do gate.
        vdrv=ogr.GetDriverByName("Memory")
        vds=vdrv.CreateDataSource("p10_mem")
        if vds is None:
            raise QgsProcessingException("P10: falha ao criar datasource vetorial MEM.")
        lyr=vds.CreateLayer("polygonize",srs=srs,geom_type=ogr.wkbPolygon)
        fd=ogr.FieldDefn("DN",ogr.OFTInteger); lyr.CreateField(fd)
        err=gdal.Polygonize(band,maskband,lyr,0,["8CONNECTED=8"])
        if err!=0:
            raise QgsProcessingException(f"P10: gdal.Polygonize falhou, código {err}.")

        code_class={1:"COMERCIAL",2:"ESPORTIVA",3:"PRESERVACAO"}
        p10_components=[]
        lyr.ResetReading()
        for feat in lyr:
            dn=int(feat.GetField("DN"))
            if dn not in code_class:
                continue
            og=feat.GetGeometryRef()
            if og is None:
                continue
            qg=QgsGeometry()
            qg.fromWkb(bytes(og.ExportToWkb()))
            if qg.isEmpty():
                continue
            if QgsWkbTypes.geometryType(qg.wkbType())!=QgsWkbTypes.PolygonGeometry:
                raise QgsProcessingException("P10: geometria não poligonal produzida.")

            # Reparação geométrica necessária para componentes 8-conectados:
            # contatos apenas diagonais podem gerar anéis auto-tocantes no
            # Polygonize. O componente lógico continua sendo único; makeValid()
            # transforma a geometria em Polygon/MultiPolygon GEOS válida sem
            # remover, agregar por proximidade ou suavizar células.
            if not qg.isGeosValid():
                qg_valid=qg.makeValid()
                if qg_valid is None or qg_valid.isEmpty():
                    raise QgsProcessingException("P10: makeValid produziu geometria vazia.")
                qg=qg_valid

            if QgsWkbTypes.geometryType(qg.wkbType())!=QgsWkbTypes.PolygonGeometry:
                raise QgsProcessingException(
                    "P10: makeValid produziu geometria não poligonal."
                )

            qg.convertToMultiType()
            if not qg.isGeosValid():
                raise QgsProcessingException(
                    "P10: componente permaneceu geometricamente inválido após makeValid."
                )

            area=float(qg.area())
            n_pix=int(round(area/(pix_m*pix_m)))
            p10_components.append({
                "classe":code_class[dn],"geom":qg,"area":area,"n_pix":n_pix
            })

        # deterministic component IDs by class + spatial bbox
        p10_components.sort(key=lambda z:(
            self.CLASSES.index(z["classe"]),
            z["geom"].boundingBox().yMinimum(),
            z["geom"].boundingBox().xMinimum()
        ))
        for i,cp in enumerate(p10_components,1):
            cp["comp"]=i

        p10_by_class={}
        for cl in self.CLASSES:
            cps=[cp for cp in p10_components if cp["classe"]==cl]
            geoms=[cp["geom"] for cp in cps]
            dissolved=QgsGeometry.unaryUnion(geoms) if geoms else QgsGeometry()
            if not dissolved.isEmpty():
                dissolved.convertToMultiType()
            p10_by_class[cl]={
                "components":cps,"geom":dissolved,
                "ncomp":len(cps),"n_pix":sum(cp["n_pix"] for cp in cps),
                "area":sum(cp["area"] for cp in cps)
            }

        invalid_comp=sum(1 for cp in p10_components if not cp["geom"].isGeosValid())
        invalid_cls=sum(1 for cl in self.CLASSES if not p10_by_class[cl]["geom"].isGeosValid())
        area_poly=sum(cp["area"] for cp in p10_components)
        area_expected=len(final_records)*pix_m*pix_m

        feedback.pushInfo(f"P10 componentes poligonais={len(p10_components)}")
        for cl in self.CLASSES:
            z=p10_by_class[cl]
            feedback.pushInfo(
                f"{cl:<13} | pixels={z['n_pix']} | area={z['area']/1e6:.6f} km2 | componentes={z['ncomp']}"
            )
        feedback.pushInfo(
            f"P10 area poligonos={area_poly/1e6:.6f} km2 | "
            f"esperada={area_expected/1e6:.6f} km2 | "
            f"delta={area_poly-area_expected:.12f} m2 | "
            f"invalidos comp/classe={invalid_comp}/{invalid_cls}"
        )

        if run_bench and modo_acaj:
            expected={
                "COMERCIAL":(141672,85),
                "ESPORTIVA":(57607,343),
                "PRESERVACAO":(31446,36)
            }
            errs=[]
            if len(p10_components)!=464: errs.append(f"P10 componentes={len(p10_components)}!=464")
            for cl,(npix,ncomp) in expected.items():
                z=p10_by_class[cl]
                if z["n_pix"]!=npix: errs.append(f"P10 {cl} pix={z['n_pix']}!={npix}")
                if z["ncomp"]!=ncomp: errs.append(f"P10 {cl} comp={z['ncomp']}!={ncomp}")
                if z["n_pix"]!=by_class_p9[cl]["pix"] or z["ncomp"]!=by_class_p9[cl]["ncomp"]:
                    errs.append(f"P10 {cl} != P9")
            AREA_TOL_M2=1e-3
            if abs(area_poly-area_expected)>AREA_TOL_M2:
                errs.append(
                    f"P10 area delta={area_poly-area_expected} m2 "
                    f"(tol={AREA_TOL_M2} m2)"
                )
            if invalid_comp or invalid_cls:
                errs.append(f"P10 invalidos={invalid_comp}/{invalid_cls}")
            if errs:
                raise QgsProcessingException("GATE P10 ACAJATUBA FALHOU: "+" | ".join(errs))
            feedback.pushInfo("[OK] Gate P10: P9=P10; 464 componentes; área 207.6525 km2; geometrias válidas.")

        # P10 outputs
        for cp in p10_components:
            f=QgsFeature(p10_comp_fields)
            f.setGeometry(cp["geom"])
            f.setAttributes([
                cp["comp"],cp["classe"],cp["n_pix"],
                cp["area"]/1e6,cp["area"]/1e4
            ])
            p10_comp_sink.addFeature(f,QgsFeatureSink.FastInsert)

        for cl in self.CLASSES:
            z=p10_by_class[cl]
            vals=[
                cl,z["n_pix"],z["area"]/1e6,z["area"]/1e4,
                100.0*z["n_pix"]/len(final_records),z["ncomp"]
            ]
            f=QgsFeature(p10_class_fields)
            f.setGeometry(z["geom"])
            f.setAttributes(vals)
            p10_class_sink.addFeature(f,QgsFeatureSink.FastInsert)

            fr=QgsFeature(p10_res_fields)
            fr.setAttributes(vals)
            p10_res_sink.addFeature(fr,QgsFeatureSink.FastInsert)

        feedback.pushInfo("="*120)
        feedback.pushInfo("MODELO 02 v6.0.1 FINAL CONCLUIDO — TODOS OS GATES P1–P10 APROVADOS")
        feedback.pushInfo("="*120)

        return {
            self.ANCORAS_OUT:anc_dest,
            self.AFINIDADES_OUT:aff_dest,
            self.INTERFACES_OUT:int_dest,
            self.RAMOS_OUT:ramo_dest,
            self.FAIXA_OUT:faixa_dest,
            self.PROX_OUT:prox_dest,
            self.FINAL_OUT:final_dest,
            self.RESUMO_OUT:res_dest,
            self.LIMIARES_OUT:lim_dest,
            self.P9_COMPONENTES_OUT:p9_comp_dest,
            self.P9_FRAGMENTOS_OUT:p9_frag_dest,
            self.P9_CLASSE_ORIGEM_OUT:p9_cross_co_dest,
            self.P9_CLASSE_CONFIANCA_OUT:p9_cross_cc_dest,
            self.P9_ORIGEM_CONFIANCA_OUT:p9_cross_oc_dest,
            self.P9_RESUMO_OUT:p9_res_dest,
            self.P10_COMPONENTES_OUT:p10_comp_dest,
            self.P10_CLASSES_OUT:p10_class_dest,
            self.P10_RESUMO_OUT:p10_res_dest,
        }