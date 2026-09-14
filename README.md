# Modelagem dos Acordos de Pesca na Amazônia

## Pipeline geoespacial para diagnóstico e segmentação de sistemas aquáticos

Este repositório apresenta um pipeline geoespacial desenvolvido para apoiar a
espacialização e a análise dos **Acordos de Pesca na Amazônia**, integrando as
coordenadas presentes nos instrumentos normativos com informações de
sensoriamento remoto sobre a distribuição e a dinâmica histórica das águas
superficiais.

O objetivo é transformar referências espaciais contidas nas Instruções
Normativas em uma representação geográfica consistente dos **sistemas aquáticos
associados às diferentes regras de manejo da pesca**, preservando as coordenadas
legais originais e incorporando informações sobre conectividade hidroespacial,
ocorrência e sazonalidade da água.

O pipeline foi desenvolvido e validado inicialmente para o **Acordo de Pesca de
Acajatuba, Amazonas, Brasil**.

> **Importante:** Acajatuba constitui o caso de desenvolvimento, calibração e
> validação do pipeline. A capacidade de generalização do método para outros
> Acordos de Pesca será avaliada por meio de aplicações independentes em outros
> sistemas aquáticos amazônicos.

---

## 1. Problema

Os instrumentos normativos que estabelecem Acordos de Pesca frequentemente
descrevem áreas e regras de manejo utilizando coordenadas geográficas associadas
a lagos, rios, paranás, canais, igarapés e outros ambientes aquáticos.

Entretanto, essas coordenadas não constituem necessariamente vértices de
polígonos e nem sempre estão posicionadas sobre o centro da lâmina d'água. Elas
podem representar margens, pontos de acesso, desembocaduras, referências locais
ou posições terrestres próximas ao ambiente regulamentado.

Por isso, a simples conexão geométrica entre coordenadas ou a aplicação de
buffers não é suficiente para representar adequadamente os sistemas aquáticos
envolvidos.

O pipeline procura responder à seguinte questão:

**Como transformar referências legais pontuais em uma representação espacial
reprodutível dos sistemas aquáticos associados às regras dos Acordos de Pesca,
considerando a configuração e a conectividade da água?**

---

## 2. Princípios metodológicos

A metodologia adota alguns princípios fundamentais:

1. **As coordenadas legais são preservadas.**  
   Os pontos extraídos das Instruções Normativas não são deslocados
   artificialmente para coincidir com a lâmina d'água.

2. **A associação ponto–água é analisada espacialmente.**  
   Distância, conectividade, ocorrência histórica e sazonalidade da água são
   utilizadas para interpretar a relação das referências legais com os sistemas
   aquáticos.

3. **Os limites funcionais não são definidos apenas por distância euclidiana.**  
   A segmentação considera afinidades hidroespaciais e zonas de competição entre
   categorias regulatórias.

4. **A automação não substitui a validação técnica.**  
   O pipeline reduz operações manuais repetitivas e explicita critérios de
   decisão, mas os resultados devem ser avaliados em relação ao instrumento
   normativo e à configuração hidrográfica local.

5. **Os produtos representam uma interpretação técnico-cartográfica.**  
   Os polígonos produzidos pelo modelo não substituem o texto legal nem
   redefinem juridicamente os limites dos Acordos de Pesca.

---

## 3. Dados de sensoriamento remoto

A principal referência utilizada para caracterização dos ambientes aquáticos é
o **JRC Global Surface Water**, produto desenvolvido pelo Joint Research Centre
da Comissão Europeia em colaboração com o Google.

O pipeline foi desenvolvido utilizando o **JRC Global Surface Water v1.4**,
derivado da série histórica Landsat e representando a dinâmica das águas
superficiais entre **1984 e 2021**, com resolução espacial de aproximadamente
**30 m**.

São utilizados três produtos:

### Maximum Water Extent — `max_extent`

Identifica os pixels nos quais foi detectada água pelo menos uma vez durante a
série histórica.

É utilizado para construir o domínio espacial potencial dos sistemas aquáticos.

### Water Occurrence — `occurrence`

Representa a frequência relativa de ocorrência de água em cada pixel.

Valores elevados indicam locais frequentemente cobertos por água, enquanto
valores menores podem representar margens, áreas sazonalmente inundadas ou
ambientes aquáticos intermitentes.

### Water Seasonality — `seasonality`

Representa o número de meses de presença de água e auxilia na diferenciação
entre ambientes permanentes e sazonais.

O uso conjunto dessas três informações permite representar melhor a
heterogeneidade hidrológica dos ambientes amazônicos do que uma imagem de uma
única data.

---

## 4. Arquitetura do pipeline

O processamento foi organizado em **dois modelos sequenciais**:

```text
INSTRUÇÕES NORMATIVAS
        │
        ▼
Extração e QA das coordenadas
        │
        ▼
Pontos legais georreferenciados
        │
        ▼
┌─────────────────────────────────────┐
│ MODELO 01                           │
│ Diagnóstico Hidroespacial           │
└─────────────────────────────────────┘
        │
        ▼
Diagnóstico JRC + domínio aquático
        │
        ▼
┌─────────────────────────────────────┐
│ MODELO 02                           │
│ Segmentação dos Sistemas Aquáticos  │
└─────────────────────────────────────┘
        │
        ▼
Afinidades e zonas de competição
        │
        ▼
Segmentação funcional
        │
        ▼
QA/QC
        │
        ▼
Polígonos dos sistemas aquáticos
        │
        ▼
VALIDAÇÃO TÉCNICA

## 5. Modelo 01 — Diagnóstico Hidroespacial

Script:
01_Diagnostico_Hidroespacial_v1.1.6.py

O Modelo 01 prepara as informações espaciais necessárias para a segmentação.

Entre suas principais operações estão:
•	reprojeção dos pontos legais para um sistema métrico adequado; 
•	geração da área operacional de processamento; 
•	recorte e preparação do JRC Maximum Water Extent; 
•	amostragem dos valores de MaxExtent, Occurrence e Seasonality; 
•	identificação dos componentes aquáticos; 
•	dissolução do domínio de água; 
•	cálculo da menor distância entre cada referência legal e o domínio aquático; 
•	geração do diagnóstico hidroespacial dos pontos; 
•	controles de qualidade das entradas e dos resultados. 
A amostragem das informações JRC é realizada preservando a grade original dos
produtos raster. A reprojeção do MaxExtent é utilizada para operações métricas
e vetoriais.

No caso de Acajatuba foi utilizado:
SIRGAS 2000 / UTM zone 20S — EPSG:31980
O CRS deve ser definido de acordo com a localização do acordo analisado.

## 6. Modelo 02 — Segmentação dos Sistemas Aquáticos

Script:
02_Segmentacao_Sistemas_Aquaticos_v2.0.5.py

O Modelo 02 utiliza os resultados do diagnóstico hidroespacial e as informações
JRC para produzir a segmentação funcional dos sistemas aquáticos.

A lógica geral envolve:

Âncoras funcionais
        ↓
Propagação hidroespacial
        ↓
Superfícies de afinidade
        ↓
Agregação por categoria
        ↓
Zonas de competição
        ↓
Interfaces entre categorias
        ↓
Ramos de equilíbrio
        ↓
Faixa operacional
        ↓
Tratamento de empates
        ↓
Segmentação funcional
        ↓
QA/QC
        ↓
Poligonização

As categorias são derivadas dos instrumentos normativos. No caso utilizado para
desenvolvimento e validação do pipeline foram analisadas categorias como:
•	COMERCIAL 
•	ESPORTIVA 
•	PRESERVAÇÃO

## 7. Afinidades e zonas de competição

O Modelo 02 não estabelece necessariamente a fronteira entre duas categorias no
ponto médio geométrico entre referências legais.

São construídas superfícies de afinidade que combinam informações sobre:
•	conectividade do domínio aquático; 
•	distância hidroespacial; 
•	características morfológicas; 
•	ocorrência histórica da água; 
•	sazonalidade; 
•	posição das referências legais. 

A comparação entre as superfícies permite identificar regiões onde duas ou mais
categorias apresentam influência semelhante.
Essas regiões constituem zonas de competição funcional.

A partir delas são identificadas interfaces e reconstruídos ramos de
equilíbrio, utilizados como elementos auxiliares na resolução das transições
entre sistemas.

Essa abordagem permite que uma fronteira seja espacialmente assimétrica quando
essa assimetria é sustentada pela configuração hidroespacial e pelas referências
normativas.

## 8. Segmentação funcional

A classificação final é realizada sobre a grade analítica de aproximadamente
30 m.

Cada pixel do domínio aquático recebe informações sobre:
•	categoria final; 
•	origem da decisão; 
•	afinidades; 
•	condição de competição; 
•	nível de confiança. 

Isso permite distinguir áreas com classificação funcional clara de regiões de
transição ou de menor confiança, que devem receber maior atenção durante a
validação técnica.

## 9. QA/QC

O pipeline inclui controles explícitos de qualidade e rastreabilidade.

Entre os controles estão:
•	conservação do domínio aquático; 
•	número de pixels processados; 
•	conectividade espacial; 
•	identificação de componentes; 
•	identificação de microfragmentos; 
•	consistência entre categoria, origem da decisão e confiança; 
•	controle das zonas de competição; 
•	validade das geometrias; 
•	fechamento entre a classificação raster e os produtos poligonais. 

A etapa de QA/QC é diagnóstica. Ela não deve modificar a classificação
funcional para melhorar artificialmente o resultado cartográfico.

## 10. Polígonos finais

Após a aprovação dos controles de QA/QC, a grade classificada é convertida em
polígonos.

A poligonização:
1.	reconstrói espacialmente as células da grade analítica; 
2.	identifica componentes conectados; 
3.	gera polígonos por componente; 
4.	dissolve os componentes por categoria; 
5.	calcula áreas e estatísticas; 
6.	verifica a consistência com a segmentação anterior. 
A poligonização é uma transformação geométrica do resultado analítico e não
realiza nova inferência sobre as categorias.

## 11. Caso de validação — Acajatuba

O Acordo de Pesca de Acajatuba, Amazonas, foi utilizado como caso de
desenvolvimento e validação metodológica.

O domínio aquático analisado apresentou:
Categoria	Pixels	Área (km²)	Área (ha)	%
COMERCIAL	141.672	127,5048	12.750,48	61,40
ESPORTIVA	57.607	51,8463	5.184,63	24,97
PRESERVAÇÃO	31.446	28,3014	2.830,14	13,63
Total	230.725	207,6525	20.765,25	100,00

A validação verificou, entre outros aspectos:
•	consistência das âncoras funcionais; 
•	propagação das afinidades; 
•	localização das zonas de competição; 
•	coerência dos ramos de equilíbrio; 
•	tratamento das zonas de transição; 
•	continuidade hidroespacial; 
•	conservação integral da classificação durante a poligonização; 
•	validade geométrica dos produtos finais. 

Acajatuba deve ser entendida como baseline de desenvolvimento e regressão,
não como demonstração suficiente de validade universal do método.

## 12. Aplicação a novos Acordos de Pesca

A próxima etapa de desenvolvimento consiste na aplicação do pipeline a outros
Acordos de Pesca da Amazônia.

Essas aplicações terão como objetivo avaliar:
•	capacidade de generalização; 
•	robustez diante de diferentes configurações hidrográficas; 
•	estabilidade das regras de segmentação; 
•	comportamento em diferentes categorias normativas; 
•	necessidade de intervenção humana; 
•	ocorrência de exceções; 
•	desempenho em áreas com diferentes graus de sazonalidade e conectividade. 

Os parâmetros não devem ser modificados apenas para produzir resultados
visualmente semelhantes aos obtidos em Acajatuba.

Resultados inesperados devem ser avaliados para distinguir:
1.	limitações do modelo; 
2.	particularidades hidrogeomorfológicas legítimas; 
3.	ambiguidades ou limitações do instrumento normativo; 
4.	problemas nos dados de entrada.

## 13. Estrutura do repositório

Modelagem-Acordos-de-Pesca-Amazonia/
│
├── README.md
│
├── scripts/
│   ├── 01_Diagnostico_Hidroespacial_v1.1.6.py
│   └── 02_Segmentacao_Sistemas_Aquaticos_v2.0.5.py
│
└── docs/
    └── documentação metodológica

Os rasters JRC, GeoPackages de processamento e demais arquivos geoespaciais de
grande volume não fazem parte do código-fonte do pipeline.

## 14. Ambiente de desenvolvimento e validação

O baseline de Acajatuba foi desenvolvida e validada no seguinte ambiente:

Componente	Versão
QGIS	3.44.8-Solothurn
Qt	5.15.13
Python	3.12.13
GDAL	3.12.2
GEOS	3.14.1
PROJ	9.8.0
PDAL	2.10.0

A execução em outras versões deve ser testada antes de ser considerada
equivalente à baseline validada.

## 15. Status

Versão atual: baseline validada em Acajatuba.
•	Extração e organização das referências legais 
•	Diagnóstico hidroespacial 
•	Integração com JRC Global Surface Water 
•	Construção das âncoras funcionais 
•	Propagação das afinidades 
•	Segmentação multiclasse 
•	Identificação das zonas de competição 
•	Reconstrução das interfaces 
•	Tratamento das zonas de transição 
•	QA/QC 
•	Poligonização 
•	Validação do caso Acajatuba 
•	Teste independente em outros Acordos de Pesca 
•	Avaliação da capacidade de generalização 
•	Consolidação de uma versão multicasos 

## 16. Limitações

O pipeline não deve ser utilizado como substituto da interpretação dos
instrumentos normativos.

A dinâmica dos sistemas aquáticos amazônicos envolve processos como:
•	migração de canais; 
•	expansão e retração sazonal da água; 
•	conexão temporária entre ambientes; 
•	inundação lateral; 
•	alterações morfológicas; 
•	diferenças entre a localização das referências legais e a configuração
contemporânea da paisagem. 

Por isso, os resultados devem ser submetidos à validação técnica,
cartográfica e normativa antes de serem utilizados em processos de gestão.

## 17. Potencial de aplicação

Uma espacialização sistemática dos Acordos de Pesca pode apoiar análises sobre:
•	extensão dos ambientes submetidos a diferentes regras de manejo; 
•	conectividade entre sistemas aquáticos; 
•	distribuição espacial das categorias de uso; 
•	monitoramento ambiental; 
•	dinâmica sazonal da água; 
•	pressão sobre recursos pesqueiros; 
•	relações entre ambientes aquáticos e uso da terra; 
•	avaliação de resultados de manejo; 
•	planejamento territorial e gestão pesqueira. 

A proposta é utilizar geoprocessamento, sensoriamento remoto e modelagem espacial
como ferramentas de apoio à interpretação e ao monitoramento dos Acordos de
Pesca, mantendo a rastreabilidade entre norma, dado espacial, processamento e
produto cartográfico.

Autoria
Antonio Oviedo
Desenvolvimento metodológico, modelagem geoespacial e implementação do pipeline.

Licença
A licença de distribuição e reutilização do código será definida após a fase de
validação independente do pipeline.
