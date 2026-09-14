# -*- coding: utf-8 -*-
"""
01_Diagnostico_Hidroespacial.py
Versão: 1.1.6
Projeto: Pipeline dos Acordos de Pesca – piloto Acajatuba
Ambiente-alvo: QGIS 3.44.x / GDAL 3.12.x / PROJ 9.8.x

BASE METODOLÓGICA
-----------------
Esta versão preserva integralmente a lógica validada do Modelo 01/V3:

Pontos legais
    -> reprojeção para CRS métrico de análise
    -> AOI operacional dissolvida (10 km por padrão)
    -> recorte do JRC MaxExtent pela AOI
    -> amostragem dos indicadores JRC preservando a grade raster original
    -> ramo vetorial/métrico do MaxExtent
    -> água vetorizada / componentes
    -> água dissolvida
    -> linhas de menor distância ponto–água
    -> diagnóstico final

PRINCÍPIOS
----------
1. As coordenadas legais NÃO são deslocadas.
2. A AOI é operacional e NÃO corresponde à área de influência.
3. MaxExtent, Occurrence e Seasonality são amostrados na grade original
   das bases JRC, evitando alteração da relação ponto–pixel.
4. A reprojeção do MaxExtent é usada apenas no ramo métrico/vetorial.
5. A associação final entre coordenada legal e ambiente aquático pertence
   ao Modelo 02 e permanece sob validação técnica humana.
6. Se a continuidade hidroespacial relevante interceptar a borda da AOI,
   a distância da AOI deve ser ampliada em nova execução.

ALTERAÇÃO DA v1.1.6 — QA/QC, SEM ALTERAÇÃO DA LÓGICA ESPACIAL
-----------------------------------------------------------------
Esta revisão NÃO altera nenhuma operação espacial do Modelo 01/V3.
Corrige exclusivamente dois pontos de QA/QC:
1. quando já existem campos JRC_MAX_*, JRC_OCC_* ou JRC_SEAS_* na camada
   de entrada, o controle passa a usar o campo de amostragem mais recente
   (adicionado ao final pela rotina native:rastersampling), e não o primeiro
   campo com o mesmo prefixo;
2. os IDs são normalizados para representação canônica, tratando 1, 1.0 e
   equivalentes como o mesmo identificador inteiro para fins de benchmark.

ALTERAÇÃO DA v1.1.5 — QA/QC, SEM ALTERAÇÃO DA LÓGICA ESPACIAL
-----------------------------------------------------------------
A lógica espacial, os parâmetros e os produtos do Modelo 01/V3 permanecem
inalterados em relação à v1.1.2. Esta revisão acrescenta somente verificações
preventivas e controles de regressão: compatibilidade estrutural das três bases
JRC, cobertura da AOI, aviso quando água vetorizada toca a borda operacional e
benchmark das distâncias-controle de Acajatuba. A v1.1.4 corrige apenas o QA/QC: usa os controles manuais consolidados (105,88; 186,87; 266,75 m), mantém o benchmark de amostragem independente do benchmark de distância e substitui o teste de borda incompatível com a API QgsGeometry por uma faixa interna de 1 pixel.

A nomenclatura das saídas permanece:

01_Pontos_Legais_Reprojetados
02_AOI_Operacional_10km
03_JRC_MaxExtent_Recortado_AOI
04_JRC_MaxExtent_Metrico
05_Pontos_Diagnostico_JRC
06_Componentes_Aquaticos_JRC
07_Agua_JRC_Dissolvida
08_Linhas_Distancia_JRC

Quando o usuário escolhe uma saída permanente, o caminho indicado pelo
usuário é respeitado. A nomenclatura padronizada é aplicada às saídas
temporárias criadas automaticamente.

HISTÓRICO v1.1.2
---------------
A v1.1.1 não reconhecia corretamente QgsProcessingOutputLayerDefinition
quando a interface enviava TEMPORARY_OUTPUT. A v1.1.2 corrige somente
essa questão de nomenclatura, sem alterar nenhuma etapa, parâmetro,
benchmark ou regra hidroespacial já validada.

BENCHMARK ACAJATUBA – V3
------------------------
JRC_MAX:
1,0,0,0,0,1,0,0,1,0,1,0,0

JRC_OCC:
94,0,0,0,0,3,0,0,69,0,1,0,0

JRC_SEAS:
12,0,0,0,0,2,0,0,4,0,0,0,0

Controles de distância consolidados:
ID 5  = 105,88 m
ID 10 = 186,87 m
ID 13 = 266,75 m
"""

from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterField,
    QgsProcessingParameterCrs,
    QgsProcessingParameterNumber,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterRasterDestination,
    QgsProcessingOutputNumber,
    QgsProcessingOutputString,
    QgsProcessingOutputLayerDefinition,
    QgsProcessingUtils,
    QgsRasterLayer,
    QgsGeometry,
    QgsCoordinateTransform,
    QgsProject,
    QgsProcessingContext,
)
import processing
import math


class DiagnosticoHidroespacial(QgsProcessingAlgorithm):

    # ------------------------------------------------------------------
    # ENTRADAS
    # ------------------------------------------------------------------

    PONTOS = "PONTOS"
    ID_FIELD = "ID_FIELD"
    CRS_ANALISE = "CRS_ANALISE"
    AOI_DIST = "AOI_DIST"
    RESOLUCAO_METRICA = "RESOLUCAO_METRICA"
    MAXEXTENT = "MAXEXTENT"
    OCCURRENCE = "OCCURRENCE"
    SEASONALITY = "SEASONALITY"
    CONECTIVIDADE_8 = "CONECTIVIDADE_8"

    # ------------------------------------------------------------------
    # SAÍDAS ESPACIAIS
    # ------------------------------------------------------------------

    PONTOS_REPROJ = "PONTOS_REPROJ"
    AOI = "AOI"
    MAX_CLIP_NATIVE = "MAX_CLIP_NATIVE"
    MAX_METRICO = "MAX_METRICO"
    PONTOS_DIAGNOSTICO = "PONTOS_DIAGNOSTICO"
    COMPONENTES_AGUA = "COMPONENTES_AGUA"
    AGUA_DISSOLVIDA = "AGUA_DISSOLVIDA"
    LINHAS_DISTANCIA = "LINHAS_DISTANCIA"

    # ------------------------------------------------------------------
    # SAÍDAS DE QA/QC
    # ------------------------------------------------------------------

    N_PONTOS_AGUA = "N_PONTOS_AGUA"
    IDS_PONTOS_AGUA = "IDS_PONTOS_AGUA"
    N_COMPONENTES = "N_COMPONENTES"
    VALIDACAO = "VALIDACAO"

    # ------------------------------------------------------------------
    # BENCHMARK ACAJATUBA V3
    # ------------------------------------------------------------------

    BENCH_MAX = {
        "1": 1,
        "2": 0,
        "3": 0,
        "4": 0,
        "5": 0,
        "6": 1,
        "7": 0,
        "8": 0,
        "9": 1,
        "10": 0,
        "11": 1,
        "12": 0,
        "13": 0,
    }

    BENCH_OCC = {
        "1": 94,
        "2": 0,
        "3": 0,
        "4": 0,
        "5": 0,
        "6": 3,
        "7": 0,
        "8": 0,
        "9": 69,
        "10": 0,
        "11": 1,
        "12": 0,
        "13": 0,
    }

    BENCH_SEAS = {
        "1": 12,
        "2": 0,
        "3": 0,
        "4": 0,
        "5": 0,
        "6": 2,
        "7": 0,
        "8": 0,
        "9": 4,
        "10": 0,
        "11": 0,
        "12": 0,
        "13": 0,
    }

    # Distâncias-controle V3 (m). Tolerância operacional de 1,0 m.
    BENCH_DIST = {
        "5": 105.88,
        "10": 186.87,
        "13": 266.75,
    }
    BENCH_DIST_TOL = 1.0

    # ------------------------------------------------------------------
    # METADADOS DO ALGORITMO
    # ------------------------------------------------------------------

    def tr(self, string):
        return QCoreApplication.translate(
            "DiagnosticoHidroespacial",
            string
        )

    def createInstance(self):
        return DiagnosticoHidroespacial()

    def name(self):
        return "01_diagnostico_hidroespacial"

    def displayName(self):
        return self.tr(
            "01 – Diagnóstico Hidroespacial v1.1.6"
        )

    def group(self):
        return self.tr(
            "Pipeline Acordos de Pesca"
        )

    def groupId(self):
        return "pipeline_acordos_pesca"

    def shortHelpString(self):
        return self.tr(
            "Modelo 01/V3 validado em Acajatuba. Reprojeta os pontos, "
            "gera AOI operacional dissolvida, recorta o MaxExtent, amostra "
            "MaxExtent/Occurrence/Seasonality preservando a grade JRC original, "
            "gera o ramo vetorial/métrico do MaxExtent e calcula linhas de "
            "menor distância dos pontos à água."
        )

    # ------------------------------------------------------------------
    # PARÂMETROS
    # ------------------------------------------------------------------

    def initAlgorithm(self, config=None):

        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.PONTOS,
                self.tr("Pontos legais"),
                [QgsProcessing.TypeVectorPoint],
            )
        )

        self.addParameter(
            QgsProcessingParameterField(
                self.ID_FIELD,
                self.tr("Campo identificador dos pontos"),
                parentLayerParameterName=self.PONTOS,
                type=QgsProcessingParameterField.Any,
            )
        )

        self.addParameter(
            QgsProcessingParameterCrs(
                self.CRS_ANALISE,
                self.tr("CRS métrico de análise"),
                defaultValue="EPSG:31980",
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.AOI_DIST,
                self.tr("Margem da AOI operacional (m)"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=10000.0,
                minValue=0.0,
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.RESOLUCAO_METRICA,
                self.tr("Resolução do MaxExtent no ramo métrico (m)"),
                type=QgsProcessingParameterNumber.Double,
                defaultValue=30.0,
                minValue=1.0,
            )
        )

        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.MAXEXTENT,
                self.tr("JRC MaxExtent – base estadual/original"),
            )
        )

        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.OCCURRENCE,
                self.tr("JRC Occurrence – base estadual/original"),
            )
        )

        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.SEASONALITY,
                self.tr("JRC Seasonality – base estadual/original"),
            )
        )

        self.addParameter(
            QgsProcessingParameterBoolean(
                self.CONECTIVIDADE_8,
                self.tr("Usar conectividade de 8 vizinhos na vetorização"),
                defaultValue=False,
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.PONTOS_REPROJ,
                self.tr("01 – Pontos legais reprojetados"),
                QgsProcessing.TypeVectorPoint,
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.AOI,
                self.tr("02 – AOI operacional"),
                QgsProcessing.TypeVectorPolygon,
            )
        )

        self.addParameter(
            QgsProcessingParameterRasterDestination(
                self.MAX_CLIP_NATIVE,
                self.tr("03 – JRC MaxExtent recortado pela AOI"),
            )
        )

        self.addParameter(
            QgsProcessingParameterRasterDestination(
                self.MAX_METRICO,
                self.tr("04 – JRC MaxExtent métrico"),
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.PONTOS_DIAGNOSTICO,
                self.tr("05 – Pontos com diagnóstico JRC"),
                QgsProcessing.TypeVectorPoint,
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.COMPONENTES_AGUA,
                self.tr("06 – Componentes aquáticos JRC"),
                QgsProcessing.TypeVectorPolygon,
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.AGUA_DISSOLVIDA,
                self.tr("07 – Água JRC dissolvida"),
                QgsProcessing.TypeVectorPolygon,
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.LINHAS_DISTANCIA,
                self.tr("08 – Linhas de distância ao JRC"),
                QgsProcessing.TypeVectorLine,
            )
        )

        self.addOutput(
            QgsProcessingOutputNumber(
                self.N_PONTOS_AGUA,
                self.tr("Número de pontos diretamente no MaxExtent")
            )
        )

        self.addOutput(
            QgsProcessingOutputString(
                self.IDS_PONTOS_AGUA,
                self.tr("IDs diretamente no MaxExtent")
            )
        )

        self.addOutput(
            QgsProcessingOutputNumber(
                self.N_COMPONENTES,
                self.tr("Número de componentes aquáticos na AOI")
            )
        )

        self.addOutput(
            QgsProcessingOutputString(
                self.VALIDACAO,
                self.tr("Resumo de QA/QC")
            )
        )

    # ------------------------------------------------------------------
    # FUNÇÕES AUXILIARES
    # ------------------------------------------------------------------

    @staticmethod
    def _sample_field(layer, prefix):
        """
        Retorna o campo de amostragem mais recente para o prefixo.

        No QGIS 3.44, native:rastersampling acrescenta o novo campo ao final
        da tabela. Se a camada de entrada já possuir, por exemplo,
        JRC_MAX_1, uma nova amostragem poderá criar JRC_MAX_2. O QA/QC não
        pode selecionar simplesmente o primeiro campo com o prefixo.
        """
        names = [
            field.name()
            for field in layer.fields()
        ]

        candidates = [
            name
            for name in names
            if name.startswith(prefix)
        ]

        if not candidates:
            raise QgsProcessingException(
                f"Campo de amostragem com prefixo "
                f"'{prefix}' não encontrado. "
                f"Campos disponíveis: {', '.join(names)}"
            )

        return candidates[-1]

    @staticmethod
    def _normalize_id(value):
        """
        Normaliza IDs apenas para comparação de QA/QC.

        Exemplos equivalentes:
        1, 1.0, '1', '1.0' -> '1'
        """
        if value is None:
            return ""

        try:
            fv = float(value)
            if math.isfinite(fv) and abs(fv - round(fv)) <= 1e-9:
                return str(int(round(fv)))
        except Exception:
            pass

        return str(value).strip()

    @staticmethod
    def _numeric_equal(a, b, tol=1e-6):
        if a is None:
            return False

        try:
            av = float(a)
            bv = float(b)

            return (
                math.isfinite(av)
                and abs(av - bv) <= tol
            )

        except Exception:
            return False

    @staticmethod
    def _named_output(parameters, key, filename):
        """
        Retorna um destino com nome legível para saídas temporárias.

        O QGIS pode encapsular TEMPORARY_OUTPUT em
        QgsProcessingOutputLayerDefinition. A v1.1.1 não tratava esse
        caso, por isso os algoritmos filhos continuavam produzindo nomes
        automáticos como Reprojetado_a_, Bordeada_, Encontrando_feicoes_,
        Dissolvido_ e Linhas_mais_curtas_.

        A v1.1.2 extrai corretamente o atributo .sink. Se a saída for
        temporária, cria um arquivo temporário com basename padronizado.
        Se o usuário definir uma saída permanente, o destino é respeitado.
        """

        value = parameters.get(key)

        if isinstance(value, QgsProcessingOutputLayerDefinition):
            sink = value.sink
        else:
            sink = value

        sink_text = "" if sink is None else str(sink)

        if (
            sink is None
            or sink_text == ""
            or sink_text.upper() == "TEMPORARY_OUTPUT"
        ):
            return QgsProcessingUtils.generateTempFilename(filename)

        return sink

    # ------------------------------------------------------------------
    # EXECUÇÃO
    # ------------------------------------------------------------------

    def processAlgorithm(
        self,
        parameters,
        context,
        feedback
    ):

        # ==============================================================
        # LEITURA DAS ENTRADAS
        # ==============================================================

        pontos = self.parameterAsSource(
            parameters,
            self.PONTOS,
            context
        )

        id_field = self.parameterAsString(
            parameters,
            self.ID_FIELD,
            context
        )

        crs_analise = self.parameterAsCrs(
            parameters,
            self.CRS_ANALISE,
            context
        )

        aoi_dist = self.parameterAsDouble(
            parameters,
            self.AOI_DIST,
            context
        )

        res_m = self.parameterAsDouble(
            parameters,
            self.RESOLUCAO_METRICA,
            context
        )

        maxextent = self.parameterAsRasterLayer(
            parameters,
            self.MAXEXTENT,
            context
        )

        occurrence = self.parameterAsRasterLayer(
            parameters,
            self.OCCURRENCE,
            context
        )

        seasonality = self.parameterAsRasterLayer(
            parameters,
            self.SEASONALITY,
            context
        )

        conn8 = self.parameterAsBool(
            parameters,
            self.CONECTIVIDADE_8,
            context
        )

        # ==============================================================
        # VALIDAÇÃO DAS ENTRADAS
        # ==============================================================

        if pontos is None:
            raise QgsProcessingException(
                "Camada de pontos legais inválida."
            )

        if not crs_analise.isValid():
            raise QgsProcessingException(
                "CRS de análise inválido."
            )

        if (
            maxextent is None
            or occurrence is None
            or seasonality is None
        ):
            raise QgsProcessingException(
                "Uma ou mais bases JRC são inválidas."
            )

        crs_max = maxextent.crs()
        crs_occ = occurrence.crs()
        crs_seas = seasonality.crs()

        if (
            not crs_max.isValid()
            or not crs_occ.isValid()
            or not crs_seas.isValid()
        ):
            raise QgsProcessingException(
                "Uma ou mais bases JRC não possuem CRS válido."
            )

        if (
            crs_max != crs_occ
            or crs_max != crs_seas
        ):
            raise QgsProcessingException(
                "MaxExtent, Occurrence e Seasonality "
                "não estão no mesmo CRS. "
                "Use as bases estaduais/originais correspondentes "
                "ao mesmo conjunto JRC."
            )

        # QA/QC estrutural: as três bases devem compartilhar a mesma
        # grade espacial original (dimensões, resolução e extensão).
        rasters = [
            ("MaxExtent", maxextent),
            ("Occurrence", occurrence),
            ("Seasonality", seasonality),
        ]

        ref_w = maxextent.width()
        ref_h = maxextent.height()
        ref_px = abs(maxextent.rasterUnitsPerPixelX())
        ref_py = abs(maxextent.rasterUnitsPerPixelY())
        ref_ext = maxextent.extent()
        tol_grid = 1e-9

        for rname, rlayer in rasters[1:]:
            px = abs(rlayer.rasterUnitsPerPixelX())
            py = abs(rlayer.rasterUnitsPerPixelY())
            ext = rlayer.extent()
            same_shape = (rlayer.width() == ref_w and rlayer.height() == ref_h)
            same_res = (abs(px-ref_px) <= tol_grid and abs(py-ref_py) <= tol_grid)
            same_extent = (
                abs(ext.xMinimum()-ref_ext.xMinimum()) <= tol_grid and
                abs(ext.xMaximum()-ref_ext.xMaximum()) <= tol_grid and
                abs(ext.yMinimum()-ref_ext.yMinimum()) <= tol_grid and
                abs(ext.yMaximum()-ref_ext.yMaximum()) <= tol_grid
            )
            if not (same_shape and same_res and same_extent):
                raise QgsProcessingException(
                    f"Grade JRC incompatível: {rname} não coincide com o "
                    "MaxExtent. Use o mesmo conjunto estadual/original."
                )

        # ==============================================================
        # DESTINOS PADRONIZADOS
        # ==============================================================

        aoi_km = aoi_dist / 1000.0

        if abs(aoi_km - round(aoi_km)) < 1e-9:
            aoi_txt = f"{int(round(aoi_km))}km"
        else:
            aoi_txt = (
                f"{aoi_km:.1f}km"
                .replace(".", "p")
            )

        out_pontos_reproj = self._named_output(
            parameters,
            self.PONTOS_REPROJ,
            "01_Pontos_Legais_Reprojetados.gpkg"
        )

        out_aoi = self._named_output(
            parameters,
            self.AOI,
            f"02_AOI_Operacional_{aoi_txt}.gpkg"
        )

        out_max_clip = self._named_output(
            parameters,
            self.MAX_CLIP_NATIVE,
            "03_JRC_MaxExtent_Recortado_AOI.tif"
        )

        out_max_metrico = self._named_output(
            parameters,
            self.MAX_METRICO,
            "04_JRC_MaxExtent_Metrico.tif"
        )

        out_pontos_diag = self._named_output(
            parameters,
            self.PONTOS_DIAGNOSTICO,
            "05_Pontos_Diagnostico_JRC.gpkg"
        )

        out_componentes = self._named_output(
            parameters,
            self.COMPONENTES_AGUA,
            "06_Componentes_Aquaticos_JRC.gpkg"
        )

        out_agua_diss = self._named_output(
            parameters,
            self.AGUA_DISSOLVIDA,
            "07_Agua_JRC_Dissolvida.gpkg"
        )

        out_linhas = self._named_output(
            parameters,
            self.LINHAS_DISTANCIA,
            "08_Linhas_Distancia_JRC.gpkg"
        )

        # ==============================================================
        # LOG INICIAL
        # ==============================================================

        feedback.pushInfo(
            "=" * 96
        )

        feedback.pushInfo(
            "DIAGNÓSTICO HIDROESPACIAL – "
            "v1.1.6 | MODELO 01/V3"
        )

        feedback.pushInfo(
            "=" * 96
        )

        feedback.pushInfo(
            f"CRS dos pontos de entrada: "
            f"{pontos.sourceCrs().authid()}"
        )

        feedback.pushInfo(
            f"CRS JRC original: "
            f"{crs_max.authid()}"
        )

        feedback.pushInfo(
            f"CRS métrico de análise: "
            f"{crs_analise.authid()}"
        )

        feedback.pushInfo(
            f"AOI operacional: "
            f"{aoi_dist:.0f} m"
        )

        feedback.pushInfo(
            f"Resolução do ramo métrico: "
            f"{res_m:.2f} m"
        )

        feedback.pushInfo(
            "Regra de amostragem: valores JRC serão lidos "
            "na grade original; a reprojeção será usada apenas "
            "no ramo métrico/vetorial."
        )

        # ==============================================================
        # 1. PONTOS NO CRS MÉTRICO
        # ==============================================================

        feedback.setProgress(5)

        feedback.pushInfo(
            "\n[1/11] Reprojetando pontos "
            "para o CRS métrico..."
        )

        pontos_reproj_result = processing.run(
            "native:reprojectlayer",
            {
                "INPUT":
                    parameters[self.PONTOS],

                "TARGET_CRS":
                    crs_analise,

                "CONVERT_CURVED_GEOMETRIES":
                    False,

                "OUTPUT":
                    out_pontos_reproj,
            },
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        pontos_reproj = (
            pontos_reproj_result["OUTPUT"]
        )

        # ==============================================================
        # 2. AOI OPERACIONAL DISSOLVIDA
        # ==============================================================

        feedback.setProgress(12)

        feedback.pushInfo(
            "\n[2/11] Gerando AOI operacional dissolvida..."
        )

        aoi_result = processing.run(
            "native:buffer",
            {
                "INPUT":
                    pontos_reproj,

                "DISTANCE":
                    aoi_dist,

                "SEGMENTS":
                    20,

                "END_CAP_STYLE":
                    0,

                "JOIN_STYLE":
                    0,

                "MITER_LIMIT":
                    2,

                "DISSOLVE":
                    True,

                "SEPARATE_DISJOINT":
                    False,

                "OUTPUT":
                    out_aoi,
            },
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        aoi = (
            aoi_result["OUTPUT"]
        )

        # QA/QC de cobertura: a AOI precisa estar integralmente coberta
        # pelas três bases JRC originais.
        aoi_layer = QgsProcessingUtils.mapLayerFromString(aoi, context)
        if aoi_layer is None:
            raise QgsProcessingException(
                "Falha ao carregar a AOI para QA/QC de cobertura."
            )

        aoi_extent = aoi_layer.extent()
        try:
            tr_to_jrc = QgsCoordinateTransform(
                crs_analise, crs_max, QgsProject.instance()
            )
            aoi_extent_jrc = tr_to_jrc.transformBoundingBox(aoi_extent)
        except Exception as exc:
            raise QgsProcessingException(
                f"Falha ao transformar a AOI para o CRS JRC: {exc}"
            )

        for rname, rlayer in rasters:
            rext = rlayer.extent()
            covered = (
                aoi_extent_jrc.xMinimum() >= rext.xMinimum() and
                aoi_extent_jrc.xMaximum() <= rext.xMaximum() and
                aoi_extent_jrc.yMinimum() >= rext.yMinimum() and
                aoi_extent_jrc.yMaximum() <= rext.yMaximum()
            )
            if not covered:
                raise QgsProcessingException(
                    f"A AOI ultrapassa a cobertura da base JRC {rname}. "
                    "Verifique os tiles/base estadual antes de continuar."
                )

        # --------------------------------------------------------------
        # AOI no CRS do raster para recorte sem alterar a grade
        # --------------------------------------------------------------

        aoi_jrc = processing.run(
            "native:reprojectlayer",
            {
                "INPUT":
                    aoi,

                "TARGET_CRS":
                    crs_max,

                "CONVERT_CURVED_GEOMETRIES":
                    False,

                "OUTPUT":
                    QgsProcessing.TEMPORARY_OUTPUT,
            },
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )["OUTPUT"]

        # ==============================================================
        # 3. RECORTE DO MAXEXTENT PELA AOI
        #    PRESERVANDO A GRADE JRC
        # ==============================================================

        feedback.setProgress(22)

        feedback.pushInfo(
            "\n[3/11] Recortando MaxExtent pela AOI "
            "sem alterar a grade original..."
        )

        max_clip_result = processing.run(
            "gdal:cliprasterbymasklayer",
            {
                "INPUT":
                    maxextent,

                "MASK":
                    aoi_jrc,

                "SOURCE_CRS":
                    None,

                "TARGET_CRS":
                    None,

                "TARGET_EXTENT":
                    None,

                "NODATA":
                    -9999,

                "ALPHA_BAND":
                    False,

                "CROP_TO_CUTLINE":
                    True,

                "KEEP_RESOLUTION":
                    True,

                "SET_RESOLUTION":
                    False,

                "X_RESOLUTION":
                    None,

                "Y_RESOLUTION":
                    None,

                "MULTITHREADING":
                    True,

                "OPTIONS":
                    "",

                "DATA_TYPE":
                    2,

                "EXTRA":
                    "-srcnodata None",

                "OUTPUT":
                    out_max_clip,
            },
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        max_clip_native = (
            max_clip_result["OUTPUT"]
        )

        # ==============================================================
        # 4. PONTOS NO CRS JRC PARA AMOSTRAGEM NATIVA
        # ==============================================================

        feedback.setProgress(30)

        feedback.pushInfo(
            "\n[4/11] Preparando pontos para amostragem "
            "na grade JRC original..."
        )

        pontos_jrc = processing.run(
            "native:reprojectlayer",
            {
                "INPUT":
                    parameters[self.PONTOS],

                "TARGET_CRS":
                    crs_max,

                "CONVERT_CURVED_GEOMETRIES":
                    False,

                "OUTPUT":
                    QgsProcessing.TEMPORARY_OUTPUT,
            },
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )["OUTPUT"]

        # ==============================================================
        # 5. AMOSTRAGEM MAXEXTENT NA BASE ORIGINAL
        # ==============================================================

        feedback.setProgress(36)

        feedback.pushInfo(
            "\n[5/11] Amostrando JRC MaxExtent "
            "na grade original..."
        )

        pontos_max = processing.run(
            "native:rastersampling",
            {
                "INPUT":
                    pontos_jrc,

                "RASTERCOPY":
                    maxextent,

                "COLUMN_PREFIX":
                    "JRC_MAX_",

                "OUTPUT":
                    QgsProcessing.TEMPORARY_OUTPUT,
            },
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )["OUTPUT"]

        # ==============================================================
        # 6. AMOSTRAGEM OCCURRENCE NA BASE ORIGINAL
        # ==============================================================

        feedback.setProgress(42)

        feedback.pushInfo(
            "\n[6/11] Amostrando JRC Occurrence "
            "na grade original..."
        )

        pontos_occ = processing.run(
            "native:rastersampling",
            {
                "INPUT":
                    pontos_max,

                "RASTERCOPY":
                    occurrence,

                "COLUMN_PREFIX":
                    "JRC_OCC_",

                "OUTPUT":
                    QgsProcessing.TEMPORARY_OUTPUT,
            },
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )["OUTPUT"]

        # ==============================================================
        # 7. AMOSTRAGEM SEASONALITY NA BASE ORIGINAL
        # ==============================================================

        feedback.setProgress(48)

        feedback.pushInfo(
            "\n[7/11] Amostrando JRC Seasonality "
            "na grade original..."
        )

        pontos_seas = processing.run(
            "native:rastersampling",
            {
                "INPUT":
                    pontos_occ,

                "RASTERCOPY":
                    seasonality,

                "COLUMN_PREFIX":
                    "JRC_SEAS_",

                "OUTPUT":
                    QgsProcessing.TEMPORARY_OUTPUT,
            },
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )["OUTPUT"]

        # --------------------------------------------------------------
        # Diagnóstico final no CRS métrico
        # mantendo os valores já amostrados
        # --------------------------------------------------------------

        pontos_diag_result = processing.run(
            "native:reprojectlayer",
            {
                "INPUT":
                    pontos_seas,

                "TARGET_CRS":
                    crs_analise,

                "CONVERT_CURVED_GEOMETRIES":
                    False,

                "OUTPUT":
                    out_pontos_diag,
            },
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        pontos_diag = (
            pontos_diag_result["OUTPUT"]
        )

        # ==============================================================
        # 8. RAMO MÉTRICO DO MAXEXTENT RECORTADO
        # ==============================================================

        feedback.setProgress(57)

        feedback.pushInfo(
            "\n[8/11] Construindo ramo métrico "
            "do MaxExtent recortado..."
        )

        max_metric_result = processing.run(
            "gdal:warpreproject",
            {
                "INPUT":
                    max_clip_native,

                "SOURCE_CRS":
                    None,

                "TARGET_CRS":
                    crs_analise,

                # nearest neighbour
                "RESAMPLING":
                    0,

                "NODATA":
                    -9999,

                "TARGET_RESOLUTION":
                    res_m,

                "OPTIONS":
                    "",

                # Int16
                "DATA_TYPE":
                    2,

                "TARGET_EXTENT":
                    None,

                "TARGET_EXTENT_CRS":
                    None,

                "MULTITHREADING":
                    True,

                "EXTRA":
                    "-tap -srcnodata -9999",

                "OUTPUT":
                    out_max_metrico,
            },
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        max_metrico = (
            max_metric_result["OUTPUT"]
        )

        max_metric_layer = QgsRasterLayer(
            max_metrico,
            "04_JRC_MaxExtent_Metrico"
        )

        if not max_metric_layer.isValid():
            raise QgsProcessingException(
                "Falha ao carregar o MaxExtent métrico."
            )

        stats = processing.run(
            "native:rasterlayerstatistics",
            {
                "INPUT":
                    max_metrico,

                "BAND":
                    1,
            },
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        feedback.pushInfo(
            "MaxExtent métrico: "
            f"MIN={float(stats['MIN']):.0f}, "
            f"MAX={float(stats['MAX']):.0f}, "
            f"MEAN={float(stats['MEAN']):.6f}"
        )

        # ==============================================================
        # 9. VETORIZAÇÃO E EXTRAÇÃO DA ÁGUA
        # ==============================================================

        feedback.setProgress(68)

        feedback.pushInfo(
            "\n[9/11] Vetorizando MaxExtent da AOI "
            "e extraindo componentes de água..."
        )

        polygonized = processing.run(
            "gdal:polygonize",
            {
                "INPUT":
                    max_metrico,

                "BAND":
                    1,

                "FIELD":
                    "JRC_MAX",

                "EIGHT_CONNECTEDNESS":
                    conn8,

                "EXTRA":
                    "",

                "OUTPUT":
                    QgsProcessing.TEMPORARY_OUTPUT,
            },
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )["OUTPUT"]

        agua_result = processing.run(
            "native:extractbyexpression",
            {
                "INPUT":
                    polygonized,

                "EXPRESSION":
                    '"JRC_MAX" >= 0.5 '
                    'AND "JRC_MAX" < 1.5',

                "OUTPUT":
                    out_componentes,
            },
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        componentes_agua = (
            agua_result["OUTPUT"]
        )

        comp_layer = (
            QgsProcessingUtils.mapLayerFromString(
                componentes_agua,
                context
            )
        )

        if comp_layer is None:
            raise QgsProcessingException(
                "Falha ao carregar os componentes aquáticos."
            )

        n_componentes = (
            comp_layer.featureCount()
        )

        feedback.pushInfo(
            f"Componentes aquáticos dentro da AOI: "
            f"{n_componentes}"
        )

        # QA/QC: água que toca a borda da AOI exige inspeção.
        try:
            aoi_geom = None
            for _f in aoi_layer.getFeatures():
                g = _f.geometry()
                if g is None or g.isEmpty():
                    continue
                aoi_geom = QgsGeometry(g) if aoi_geom is None else aoi_geom.combine(g)

            n_borda = 0
            if aoi_geom is not None and not aoi_geom.isEmpty():
                # QA robusto: considera a faixa interna de 1 pixel (res_m)
                # junto à borda da AOI. Não altera qualquer geometria de análise.
                inner = aoi_geom.buffer(-float(res_m), 8)
                if inner is not None and not inner.isEmpty():
                    edge_ring = aoi_geom.difference(inner)
                else:
                    edge_ring = QgsGeometry(aoi_geom)

                for _cf in comp_layer.getFeatures():
                    cg = _cf.geometry()
                    if (
                        cg is not None
                        and not cg.isEmpty()
                        and cg.intersects(edge_ring)
                    ):
                        n_borda += 1

            if n_borda > 0:
                feedback.pushWarning(
                    f"[QA] {n_borda} componente(s) aquático(s) tocam a borda "
                    "da AOI. Se houver continuidade hidroespacial relevante, "
                    "amplie a AOI e execute novamente o Modelo 01."
                )
            else:
                feedback.pushInfo(
                    "[OK] Nenhum componente aquático toca a borda da AOI."
                )
        except Exception as exc:
            feedback.pushWarning(
                f"[QA] Teste água–borda não concluído: {exc}"
            )

        # ==============================================================
        # 10. DISSOLVER ÁGUA
        # ==============================================================

        feedback.setProgress(78)

        feedback.pushInfo(
            "\n[10/11] Dissolvendo a água "
            "do MaxExtent..."
        )

        agua_diss_result = processing.run(
            "native:dissolve",
            {
                "INPUT":
                    componentes_agua,

                "FIELD":
                    [],

                "SEPARATE_DISJOINT":
                    False,

                "OUTPUT":
                    out_agua_diss,
            },
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        agua_dissolvida = (
            agua_diss_result["OUTPUT"]
        )

        # ==============================================================
        # 11. MENOR DISTÂNCIA DOS PONTOS À ÁGUA
        # ==============================================================

        feedback.setProgress(88)

        feedback.pushInfo(
            "\n[11/11] Calculando menor distância "
            "dos pontos legais à água..."
        )

        linhas_result = processing.run(
            "native:shortestline",
            {
                "SOURCE":
                    pontos_diag,

                "DESTINATION":
                    agua_dissolvida,

                "METHOD":
                    0,

                "NEIGHBORS":
                    1,

                "DISTANCE":
                    None,

                "OUTPUT":
                    out_linhas,
            },
            context=context,
            feedback=feedback,
            is_child_algorithm=True,
        )

        linhas_distancia = (
            linhas_result["OUTPUT"]
        )

        # QA/QC – DISTÂNCIAS-CONTROLE ACAJATUBA V3
        linhas_layer = QgsProcessingUtils.mapLayerFromString(
            linhas_distancia, context
        )
        dist_bench_applicable = False
        dist_bench_ok = True
        dist_rows = {}

        if linhas_layer is not None and id_field in {f.name() for f in linhas_layer.fields()}:
            ids_line = set()
            for lf in linhas_layer.getFeatures():
                pid_l = self._normalize_id(lf[id_field])
                ids_line.add(pid_l)
                if pid_l in self.BENCH_DIST:
                    geom_l = lf.geometry()
                    if geom_l is not None and not geom_l.isEmpty():
                        dist_rows[pid_l] = float(geom_l.length())

            dist_bench_applicable = set(self.BENCH_DIST).issubset(ids_line)
            if dist_bench_applicable:
                for pid_l, expected in self.BENCH_DIST.items():
                    got = dist_rows.get(pid_l)
                    if got is None or abs(got-expected) > self.BENCH_DIST_TOL:
                        dist_bench_ok = False

                if dist_bench_ok:
                    feedback.pushInfo(
                        "[OK] Distâncias-controle Acajatuba V3 reproduzidas "
                        f"(±{self.BENCH_DIST_TOL:.1f} m)."
                    )
                else:
                    feedback.pushWarning(
                        "[ATENÇÃO] Distâncias-controle Acajatuba V3 divergentes."
                    )

                for pid_l in ("5", "10", "13"):
                    if pid_l in dist_rows:
                        feedback.pushInfo(
                            f"ID {pid_l}: distância={dist_rows[pid_l]:.2f} m | "
                            f"esperado≈{self.BENCH_DIST[pid_l]:.2f} m"
                        )

        # ==============================================================
        # QA/QC – BENCHMARK ACAJATUBA
        # ==============================================================

        diag_layer = (
            QgsProcessingUtils.mapLayerFromString(
                pontos_diag,
                context
            )
        )

        if diag_layer is None:
            raise QgsProcessingException(
                "Falha ao carregar os pontos de diagnóstico."
            )

        f_max = self._sample_field(
            diag_layer,
            "JRC_MAX_"
        )

        f_occ = self._sample_field(
            diag_layer,
            "JRC_OCC_"
        )

        f_seas = self._sample_field(
            diag_layer,
            "JRC_SEAS_"
        )

        feedback.pushInfo(
            "Campos JRC usados no QA/QC: "
            f"MAX={f_max} | OCC={f_occ} | SEAS={f_seas}"
        )

        hits = []
        benchmark_applicable = True
        benchmark_ok = True
        rows = {}

        for feat in diag_layer.getFeatures():

            pid = self._normalize_id(
                feat[id_field]
            )

            v_max = feat[f_max]
            v_occ = feat[f_occ]
            v_seas = feat[f_seas]

            rows[pid] = (
                v_max,
                v_occ,
                v_seas
            )

            if self._numeric_equal(
                v_max,
                1
            ):
                hits.append(
                    pid
                )

            if pid not in self.BENCH_MAX:
                benchmark_applicable = False

            else:

                if not self._numeric_equal(
                    v_max,
                    self.BENCH_MAX[pid]
                ):
                    benchmark_ok = False

                if not self._numeric_equal(
                    v_occ,
                    self.BENCH_OCC[pid]
                ):
                    benchmark_ok = False

                if not self._numeric_equal(
                    v_seas,
                    self.BENCH_SEAS[pid]
                ):
                    benchmark_ok = False

        try:
            hits = [
                str(value)
                for value in sorted(
                    int(item)
                    for item in hits
                )
            ]

        except Exception:
            hits = sorted(
                hits
            )

        ids_hits = ", ".join(
            hits
        )

        n_hits = len(
            hits
        )

        feedback.pushInfo(
            "-" * 96
        )

        feedback.pushInfo(
            "QA/QC – AMOSTRAGEM NA GRADE JRC ORIGINAL"
        )

        feedback.pushInfo(
            "-" * 96
        )

        feedback.pushInfo(
            f"Pontos diretamente no MaxExtent: "
            f"{n_hits}/{pontos.featureCount()}"
        )

        feedback.pushInfo(
            f"IDs: "
            f"{ids_hits if ids_hits else 'nenhum'}"
        )

        ids_presentes = set(
            rows.keys()
        )

        benchmark_applicable = (
            benchmark_applicable
            and ids_presentes
            == set(
                self.BENCH_MAX.keys()
            )
        )

        if benchmark_applicable:

            # MAX/OCC/SEAS e distâncias são controles independentes.
            # Uma divergência de distância não deve mascarar a validação
            # da amostragem na grade JRC original.
            if benchmark_ok:

                validacao = (
                    "APROVADO – amostragem reproduz "
                    "o benchmark Acajatuba V3 "
                    "(MAX/OCC/SEAS na grade original)."
                )

                feedback.pushInfo(
                    "[OK] Benchmark Acajatuba V3 reproduzido."
                )

            else:

                validacao = (
                    "VERIFICAR – amostragem não reproduziu "
                    "integralmente o benchmark Acajatuba V3. "
                    "Conferir se foram selecionadas as bases "
                    "estaduais/originais corretas."
                )

                feedback.pushInfo(
                    "[ATENÇÃO] Benchmark Acajatuba V3 divergente."
                )

                def sort_key(value):
                    if value.isdigit():
                        return int(value)
                    return value

                for pid in sorted(
                    rows,
                    key=sort_key
                ):

                    feedback.pushInfo(
                        f"ID {pid}: "
                        f"MAX={rows[pid][0]} | "
                        f"OCC={rows[pid][1]} | "
                        f"SEAS={rows[pid][2]}"
                    )

        else:

            validacao = (
                "EXECUÇÃO OK – benchmark Acajatuba "
                "não aplicado porque o conjunto de IDs "
                "não corresponde exatamente a 1–13."
            )

            feedback.pushInfo(
                "[INFO] Benchmark Acajatuba "
                "não aplicável a este conjunto de pontos."
            )

        feedback.pushInfo(
            "Observação: a contagem de componentes "
            "é diagnóstica e depende da AOI, resolução "
            "e regra de conectividade; "
            "não é usada como benchmark rígido."
        )

        feedback.pushInfo(
            "Saídas temporárias padronizadas:"
        )

        feedback.pushInfo(
            "01_Pontos_Legais_Reprojetados"
        )

        feedback.pushInfo(
            f"02_AOI_Operacional_{aoi_txt}"
        )

        feedback.pushInfo(
            "03_JRC_MaxExtent_Recortado_AOI"
        )

        feedback.pushInfo(
            "04_JRC_MaxExtent_Metrico"
        )

        feedback.pushInfo(
            "05_Pontos_Diagnostico_JRC"
        )

        feedback.pushInfo(
            "06_Componentes_Aquaticos_JRC"
        )

        feedback.pushInfo(
            "07_Agua_JRC_Dissolvida"
        )

        feedback.pushInfo(
            "08_Linhas_Distancia_JRC"
        )

        # =============================================================
        # REGISTRO EXPLÍCITO DAS OITO SAÍDAS NO PAINEL DO QGIS
        # =============================================================
        # As versões anteriores executavam corretamente os algoritmos filhos,
        # porém os destinos temporários eram retornados como caminhos/IDs internos
        # sem registro explícito para carregamento ao final do algoritmo pai.
        # Aqui todas as oito saídas são registradas nominalmente no contexto.
        # Esta correção NÃO altera nenhuma operação espacial.
        _outputs_to_load = [
            (pontos_reproj, "01_Pontos_Legais_Reprojetados"),
            (aoi, f"02_AOI_Operacional_{aoi_txt}"),
            (max_clip_native, "03_JRC_MaxExtent_Recortado_AOI"),
            (max_metrico, "04_JRC_MaxExtent_Metrico"),
            (pontos_diag, "05_Pontos_Diagnostico_JRC"),
            (componentes_agua, "06_Componentes_Aquaticos_JRC"),
            (agua_dissolvida, "07_Agua_JRC_Dissolvida"),
            (linhas_distancia, "08_Linhas_Distancia_JRC"),
        ]

        for _source, _name in _outputs_to_load:
            try:
                _details = QgsProcessingContext.LayerDetails(
                    _name,
                    context.project(),
                    _name
                )
                context.addLayerToLoadOnCompletion(
                    str(_source),
                    _details
                )
            except Exception as _exc:
                raise QgsProcessingException(
                    f"Falha ao registrar saída '{_name}' para carregamento: {_exc}"
                )

        feedback.pushInfo(
            "[OK] As oito saídas foram registradas explicitamente para "
            "carregamento no painel do QGIS."
        )

        feedback.setProgress(
            100
        )

        return {
            self.PONTOS_REPROJ:
                pontos_reproj,

            self.AOI:
                aoi,

            self.MAX_CLIP_NATIVE:
                max_clip_native,

            self.MAX_METRICO:
                max_metrico,

            self.PONTOS_DIAGNOSTICO:
                pontos_diag,

            self.COMPONENTES_AGUA:
                componentes_agua,

            self.AGUA_DISSOLVIDA:
                agua_dissolvida,

            self.LINHAS_DISTANCIA:
                linhas_distancia,

            self.N_PONTOS_AGUA:
                n_hits,

            self.IDS_PONTOS_AGUA:
                ids_hits,

            self.N_COMPONENTES:
                n_componentes,

            self.VALIDACAO:
                validacao,
        }