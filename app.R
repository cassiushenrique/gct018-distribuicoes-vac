# ================================================================
# DISTRIBUIÇÕES VAC | GCT018 - UFLA / ICTIN
# Aulas 27 a 32. Execute: shiny::runApp("app_distribuicoes_VAC.R")
# Dependências: install.packages(c("shiny", "ggplot2"))
# ================================================================

library(shiny)
library(ggplot2)

azul <- "#123063"
realce <- "#2782AD"
fmt <- function(x, dig = 5) {
  if (is.na(x)) return("—")
  if (is.infinite(x)) return(if (x > 0) "+∞" else "−∞")
  formatC(x, digits = dig, format = "fg", decimal.mark = ",", big.mark = ".")
}
math <- function(z) HTML(paste0("\\(", htmltools::htmlEscape(z), "\\)"))

# A parametrização Weibull dos slides é f(x)=gamma/beta*x^(gamma-1)*exp(-x^gamma/beta).
# Na função dweibull do R, o argumento scale deve ser beta^(1/gamma), NÃO beta.
# A Gama dos slides é oferecida nas duas parametrizações: escala beta e taxa lambda.
criar_modelo <- function(tipo, v) {
  switch(tipo,
    uniforme = {
      a <- v$a; b <- v$b
      list(nome = "Uniforme contínua", notacao = "X\\sim U(a,b)",
           formula = "f(x)=\\begin{cases}\\dfrac{1}{b-a},&a\\leq x\\leq b,\\\\0,&\\text{caso contrário},\\end{cases}",
           dominio = "a,b\\in\\mathbb{R},\\ a<b;\\quad x\\in[a,b]",
           parametros = "a: limite inferior; b: limite superior.",
           aplicacao = "Representa valores igualmente densos em um intervalo conhecido, por exemplo, um instante de chegada aleatório em uma janela de tempo.",
           dens = function(x) dunif(x, a, b), acum = function(x) punif(x, a, b),
           quant = function(p) qunif(p, a, b), media = (a+b)/2,
           variancia = (b-a)^2/12, suporte = c(a,b), observacao = "A densidade é constante no intervalo, mas P(X = x) = 0 para qualquer ponto isolado.")
    },
    normal = {
      mu <- v$mu; sig <- v$sigma
      list(nome = "Normal", notacao = "X\\sim N(\\mu,\\sigma^2)",
           formula = "f(x)=\\dfrac{1}{\\sigma\\sqrt{2\\pi}}\\exp\\!\\left[-\\dfrac{(x-\\mu)^2}{2\\sigma^2}\\right]",
           dominio = "\\mu\\in\\mathbb{R},\\ \\sigma>0;\\quad x\\in\\mathbb{R}",
           parametros = "μ: média (localização); σ: desvio-padrão (dispersão); σ²: variância.",
           aplicacao = "Modela, sob pressupostos adequados, erros de medição e características dimensionais de processos industriais. A padronização é Z = (X − μ)/σ.",
           dens = function(x) dnorm(x, mu, sig), acum = function(x) pnorm(x, mu, sig),
           quant = function(p) qnorm(p, mu, sig), media = mu,
           variancia = sig^2, suporte = c(-Inf, Inf),
           observacao = "Alterar μ desloca a curva; alterar σ altera sua dispersão e altura, preservando área total igual a 1.")
    },
    gama_escala = {
      al <- v$alpha; be <- v$beta
      list(nome = "Gama — forma e escala", notacao = "X\\sim\\Gamma(\\alpha,\\beta)",
           formula = "f(x)=\\dfrac{x^{\\alpha-1}e^{-x/\\beta}}{\\Gamma(\\alpha)\\beta^{\\alpha}}",
           dominio = "\\alpha>0,\\ \\beta>0;\\quad x>0",
           parametros = "α: forma; β: escala (β = 1/λ).",
           aplicacao = "Modela grandezas positivas e tempos até eventos em confiabilidade. Para α inteiro, obtém-se a família Erlang; para α = 1, a exponencial.",
           dens = function(x) dgamma(x, shape=al, scale=be),
           acum = function(x) pgamma(x, shape=al, scale=be),
           quant = function(p) qgamma(p, shape=al, scale=be),
           media = al*be, variancia = al*be^2, suporte = c(0,Inf),
           observacao = "Nesta parametrização: E(X) = αβ e Var(X) = αβ².")
    },
    gama_taxa = {
      al <- v$alpha; la <- v$lambda
      list(nome = "Gama — forma e taxa", notacao = "X\\sim\\Gamma(\\alpha,\\lambda)\\;\\;[\\text{taxa}]",
           formula = "f(x)=\\dfrac{\\lambda^{\\alpha}}{\\Gamma(\\alpha)}x^{\\alpha-1}e^{-\\lambda x}",
           dominio = "\\alpha>0,\\ \\lambda>0;\\quad x>0",
           parametros = "α: forma; λ: taxa (λ = 1/β).",
           aplicacao = "Alternativa à parametrização por escala para tempos de espera e dados positivos; permite explorar o efeito de uma taxa maior sobre a concentração dos tempos.",
           dens = function(x) dgamma(x, shape=al, rate=la),
           acum = function(x) pgamma(x, shape=al, rate=la),
           quant = function(p) qgamma(p, shape=al, rate=la),
           media = al/la, variancia = al/la^2, suporte = c(0,Inf),
           observacao = "Cuidado: aqui E(X) = α/λ. Taxa λ e escala β são recíprocas, não intercambiáveis.")
    },
    beta = {
      al <- v$alpha; be <- v$beta
      list(nome = "Beta", notacao = "X\\sim\\operatorname{Beta}(\\alpha,\\beta)",
           formula = "f(x)=\\dfrac{x^{\\alpha-1}(1-x)^{\\beta-1}}{B(\\alpha,\\beta)}",
           dominio = "\\alpha>0,\\ \\beta>0;\\quad 0<x<1",
           parametros = "α e β: parâmetros de forma; B(α,β) = Γ(α)Γ(β)/Γ(α+β).",
           aplicacao = "Adequada para modelar proporções e frações contínuas estritamente entre 0 e 1, desde que não existam massas positivas nos extremos.",
           dens = function(x) dbeta(x, al, be), acum = function(x) pbeta(x, al, be),
           quant = function(p) qbeta(p, al, be),
           media = al/(al+be), variancia = al*be/((al+be)^2*(al+be+1)), suporte = c(0,1),
           observacao = "Quando α ou β é inferior a 1, a densidade pode divergir em 0 ou 1. Isso não invalida o modelo: a área total continua sendo 1.")
    },
    exponencial = {
      be <- v$beta
      list(nome = "Exponencial — escala", notacao = "X\\sim\\operatorname{Exp}(\\beta)",
           formula = "f(x)=\\dfrac{1}{\\beta}e^{-x/\\beta}",
           dominio = "\\beta>0;\\quad x>0",
           parametros = "β: escala e média; taxa de ocorrência λ = 1/β.",
           aplicacao = "Modelo de tempo até falha sob taxa de falha constante. É o caso particular da Gama com α = 1 e possui propriedade de ausência de memória.",
           dens = function(x) dexp(x, rate=1/be), acum = function(x) pexp(x, rate=1/be),
           quant = function(p) qexp(p, rate=1/be), media = be,
           variancia = be^2, suporte = c(0,Inf),
           observacao = "P(X > s+t | X > s) = P(X > t), para s,t ≥ 0.")
    },
    weibull = {
      ga <- v$gamma; be <- v$beta; escalaR <- be^(1/ga)
      list(nome = "Weibull — parametrização dos slides", notacao = "X\\sim\\operatorname{Weibull}(\\gamma,\\beta)",
           formula = "f(x)=\\dfrac{\\gamma}{\\beta}x^{\\gamma-1}\\exp\\!\\left(-\\dfrac{x^{\\gamma}}{\\beta}\\right)",
           dominio = "\\gamma>0,\\ \\beta>0;\\quad x>0",
           parametros = "γ: forma; β: parâmetro positivo na expressão exp(−x^γ/β); escala convencional do R = β^(1/γ).",
           aplicacao = "Modela tempos até falha de componentes: γ < 1 representa taxa de falha decrescente, γ = 1 constante e γ > 1 crescente.",
           dens = function(x) dweibull(x, shape=ga, scale=escalaR),
           acum = function(x) pweibull(x, shape=ga, scale=escalaR),
           quant = function(p) qweibull(p, shape=ga, scale=escalaR),
           media = escalaR*gamma(1+1/ga),
           variancia = escalaR^2*(gamma(1+2/ga)-gamma(1+1/ga)^2),
           suporte = c(0,Inf),
           observacao = "Parametrização fiel ao slide: β NÃO é o argumento scale de dweibull(); scale = β^(1/γ).")
    },
    quiquadrado = {
      gl <- v$gl
      list(nome = "Qui-quadrado", notacao = "X\\sim\\chi^2_p",
           formula = "f(x)=\\dfrac{x^{p/2-1}e^{-x/2}}{2^{p/2}\\Gamma(p/2)}",
           dominio = "p\\in\\{1,2,3,\\ldots\\};\\quad x>0",
           parametros = "p: graus de liberdade (inteiro positivo); caso particular da Gama com α = p/2 e β = 2.",
           aplicacao = "Surge na inferência da variância populacional normal e em estatísticas de testes de hipóteses.",
           dens = function(x) dchisq(x, df=gl), acum = function(x) pchisq(x, df=gl),
           quant = function(p) qchisq(p, df=gl), media = gl,
           variancia = 2*gl, suporte = c(0,Inf),
           observacao = "Distribuição incluída por aparecer explicitamente como caso particular da Gama nos slides.")
    }
  )
}

parametros_ui <- function(tipo) {
  slider <- function(id, rotulo, mn, mx, val, passo) {
    sliderInput(id, rotulo, min=mn, max=mx, value=val, step=passo,
                sep="", ticks=FALSE)
  }
  switch(tipo,
    uniforme = tagList(slider("a", "a — limite inferior", -20, 19, 0, 0.5),
                       slider("b", "b — limite superior", -19, 20, 5, 0.5),
                       helpText("Restrição: a < b. Ajuste os limites para manter a ordem.")),
    normal = tagList(slider("mu", "μ — média", -50, 50, 0, 0.5),
                     slider("sigma", "σ — desvio-padrão", 0.1, 30, 1, 0.1)),
    gama_escala = tagList(slider("alpha", "α — forma", 0.2, 15, 2, 0.1),
                          slider("beta", "β — escala", 0.1, 20, 2, 0.1)),
    gama_taxa = tagList(slider("alpha", "α — forma", 0.2, 15, 2, 0.1),
                        slider("lambda", "λ — taxa (= 1/β)", 0.1, 15, 1, 0.1)),
    beta = tagList(slider("alpha", "α — forma", 0.2, 15, 2, 0.1),
                   slider("beta", "β — forma", 0.2, 15, 5, 0.1)),
    exponencial = slider("beta", "β — escala e média", 0.1, 30, 2, 0.1),
    weibull = tagList(slider("gamma", "γ — forma", 0.2, 8, 1.5, 0.1),
                      slider("beta", "β — parâmetro dos slides", 0.1, 30, 2, 0.1)),
    quiquadrado = slider("gl", "p — graus de liberdade", 1, 30, 4, 1)
  )
}

ui <- fluidPage(
  tags$head(
    tags$meta(name="viewport", content="width=device-width, initial-scale=1"),
    tags$style(HTML("\n      body {font-family:'Segoe UI',Arial,sans-serif; background:#F3F6FB; color:#192638;}\n      .container-fluid {max-width:1560px; padding:22px 30px;}\n      .top {background:linear-gradient(105deg,#102951,#17477A); color:white;\n            border-radius:14px; padding:22px 28px; margin-bottom:20px;}\n      .top h1 {font-weight:750; font-size:31px; margin:0 0 6px;}\n      .top p {margin:0; opacity:.88;}\n      .panel-vac {background:white; border:1px solid #DEE6F0; border-radius:12px;\n                  padding:22px; margin-bottom:18px; box-shadow:0 3px 14px #142F5210;}\n      .panel-vac h3 {font-size:19px; color:#123063; margin-top:0; font-weight:700;}\n      .form-group {margin-bottom:20px;} label {color:#263C57; font-weight:600;}\n      .irs--shiny .irs-bar,.irs--shiny .irs-single {background:#2782AD; border-color:#2782AD;}\n      .irs--shiny .irs-handle {background:white; border:2px solid #2782AD;}\n      .stat {background:#F5F8FC; border:1px solid #E0E8F2; border-radius:10px;\n             padding:13px 16px; margin-bottom:12px;}\n      .stat .k {font-size:12px; text-transform:uppercase; letter-spacing:.06em; color:#59708A;}\n      .stat .v {font-size:22px; font-weight:750; color:#123063; margin-top:4px;}\n      .answer {background:#EAF5FB; border-left:5px solid #2782AD; padding:16px;\n               border-radius:8px; font-size:19px; color:#123063; margin-top:14px;}\n      .model-text {font-size:15px; line-height:1.7;}\n      .note {font-size:13px; color:#526780; line-height:1.55;}\n      .footer {font-size:12px; color:#60718A; padding:8px 3px;}\n      .btn-primary {background:#123063; border-color:#123063;}\n      @media(max-width:900px){.container-fluid{padding:12px;} .top h1{font-size:25px;}}\n    "))
  ),
  div(class="top", h1("Distribuições VAC"),
      p("GCT018 · Estatística: Conceitos e Experimentações  |  UFLA · ICTIN · Aulas 27 a 32")),
  fluidRow(
    column(4,
      div(class="panel-vac",
          h3("1. Modelo probabilístico"),
          selectInput("tipo", "Distribuição", choices=c(
            "Uniforme contínua"="uniforme", "Normal"="normal",
            "Gama — forma e escala (α, β)"="gama_escala",
            "Gama — forma e taxa (α, λ)"="gama_taxa",
            "Beta"="beta", "Exponencial"="exponencial",
            "Weibull (parametrização dos slides)"="weibull",
            "Qui-quadrado"="quiquadrado"), selected="normal"),
          uiOutput("parametros"),
          tags$hr(), h3("2. Probabilidade / quantil"),
          selectInput("evento", "Região de interesse", choices=c(
            "Cauda inferior: P(X ≤ a)"="inferior",
            "Cauda superior: P(X ≥ a)"="superior",
            "Intervalo: P(a ≤ X ≤ b)"="intervalo",
            "Região externa: P(X ≤ a ou X ≥ b)"="externa",
            "Quantil: encontrar q com P(X ≤ q) = u"="quantil"),
            selected="inferior"),
          uiOutput("limites"),
          checkboxInput("mostrar_media", "Mostrar linha da média", TRUE),
          checkboxInput("mostrar_legenda", "Mostrar limites do evento", TRUE),
          helpText("A região azul no gráfico representa a probabilidade do evento.")),
      div(class="panel-vac",
          h3("Explore"),
          p("Mude um parâmetro por vez e compare localização, dispersão, assimetria e altura da densidade."),
          p("Na distribuição contínua, probabilidades são áreas: não confunda f(x) com P(X = x)."),
          p(class="note", "Os limites do evento são redefinidos para quantis da distribuição quando os parâmetros são modificados."))
    ),
    column(8,
      div(class="panel-vac",
          h3(textOutput("titulo", inline=TRUE)),
          plotOutput("grafico", height="440px"),
          div(class="answer", uiOutput("resultado")),
          br(), fluidRow(
            column(4, div(class="stat", div(class="k", "Esperança · E(X)"), div(class="v", textOutput("media")))),
            column(4, div(class="stat", div(class="k", "Variância · Var(X)"), div(class="v", textOutput("variancia")))),
            column(4, div(class="stat", div(class="k", "Desvio-padrão · σ"), div(class="v", textOutput("dp"))))
          ),
          p(class="note", textOutput("nota_grafico"))),
      div(class="panel-vac", h3("Modelo probabilístico · definição e interpretação"),
          uiOutput("explicacao")),
      div(class="footer", "Material didático · Prof. Dr. Cássius Henrique Xavier Oliveira · UFLA / ICTIN. Base: Casella e Berger (2010) e Montgomery (2005).")
    )
  )
)

server <- function(input, output, session) {
  output$parametros <- renderUI(parametros_ui(input$tipo))

  modelo <- reactive({
    req(input$tipo)
    ids <- switch(input$tipo, uniforme=c("a","b"), normal=c("mu","sigma"),
                  gama_escala=c("alpha","beta"), gama_taxa=c("alpha","lambda"),
                  beta=c("alpha","beta"), exponencial="beta",
                  weibull=c("gamma","beta"), quiquadrado="gl")
    v <- lapply(ids, function(id) input[[id]])
    req(all(vapply(v, function(z) !is.null(z) && is.finite(z), logical(1))))
    names(v) <- ids
    validate(need(input$tipo != "uniforme" || v$a < v$b,
                  "Na uniforme, o limite inferior a deve ser estritamente menor que b."))
    criar_modelo(input$tipo, v)
  })

  # O painel de limites é reativo à distribuição e aos parâmetros.
  # Os limites seguem a escala de X (não a escala padronizada Z).
  janela <- reactive({
    m <- modelo()
    if (input$tipo == "uniforme") return(m$suporte)
    if (input$tipo == "beta") return(c(0,1))
    qs <- m$quant(c(.0005, .9995))
    if (is.finite(m$suporte[1])) qs[1] <- m$suporte[1]
    qs
  })

  output$limites <- renderUI({
    j <- janela(); m <- modelo(); req(input$evento)
    if (input$evento == "quantil")
      return(sliderInput("u", "u — probabilidade acumulada", min=.001, max=.999,
                         value=.95, step=.001, ticks=FALSE))
    passo <- max((j[2]-j[1])/1000, .000001)
    if (input$evento %in% c("intervalo", "externa")) {
      return(sliderInput("ab", "a e b — limites do evento", min=j[1], max=j[2],
                         value=as.numeric(m$quant(c(.25,.75))), step=passo,
                         ticks=FALSE, sep=""))
    }
    sliderInput("limite", "a — ponto de corte", min=j[1], max=j[2],
                value=as.numeric(m$quant(.5)), step=passo, ticks=FALSE, sep="")
  })

  evento <- reactive({
    m <- modelo(); req(input$evento)
    tipo <- input$evento
    if (tipo == "quantil") {
      req(input$u); q <- m$quant(input$u)
      return(list(a=q, b=Inf, prob=input$u, q=q, tipo=tipo,
                  descricao=paste0("q = F⁻¹(", fmt(input$u), ") = ",fmt(q),
                                   ";  P(X ≤ q) = ",fmt(input$u))))
    }
    if (tipo %in% c("intervalo", "externa")) {
      req(input$ab); a <- input$ab[1]; b <- input$ab[2]
      pa <- m$acum(a); pb <- m$acum(b)
      pr <- if (tipo == "intervalo") pb-pa else pa + m$acum(b)*0 + (1-pb)
      texto <- if (tipo == "intervalo") "P(a ≤ X ≤ b)" else "P(X ≤ a ou X ≥ b)"
      return(list(a=a,b=b,prob=pmax(0,pmin(1,pr)),tipo=tipo,
                  descricao=paste0(texto," = ",fmt(pr), "  (",fmt(100*pr),"%)")))
    }
    req(input$limite); a <- input$limite
    pr <- if (tipo == "inferior") m$acum(a) else 1-m$acum(a)
    texto <- if (tipo == "inferior") "P(X ≤ a)" else "P(X ≥ a)"
    list(a=a,b=Inf,prob=pmax(0,pmin(1,pr)),tipo=tipo,
         descricao=paste0(texto," = ",fmt(pr), "  (",fmt(100*pr),"%)"))
  })

  output$titulo <- renderText(paste("Distribuição", modelo()$nome))
  output$media <- renderText(fmt(modelo()$media, 6))
  output$variancia <- renderText(fmt(modelo()$variancia, 6))
  output$dp <- renderText(fmt(sqrt(max(0, modelo()$variancia)), 6))
  output$resultado <- renderUI({
    e <- evento()
    if (e$tipo == "quantil")
      return(HTML(paste0("<strong>Quantil:</strong> ", e$descricao)))
    HTML(paste0("<strong>Área selecionada:</strong> ",e$descricao,
                "<br><small>Limites: a = ",fmt(e$a),
                if (e$tipo %in% c("intervalo","externa")) paste0("; b = ",fmt(e$b)) else "",
                "</small>"))
  })

  output$grafico <- renderPlot({
    m <- modelo(); j <- janela(); e <- evento()
    l <- j[1]; u <- j[2]; delta <- u-l
    validate(need(is.finite(delta) && delta>0,"Limites de visualização inválidos."))
    # A borda é retirada quando uma densidade tem singularidade integrável.
    n <- 1700L
    xs <- seq(l,u,length.out=n)
    if (input$tipo == "beta") {
      xs[1] <- max(xs[1],1e-7); xs[n] <- min(xs[n],1-1e-7)
    } else if (l == 0 && input$tipo %in% c("gama_escala","gama_taxa","weibull","quiquadrado")) {
      xs[1] <- max(delta*1e-6, .Machine$double.eps)
    }
    ys <- m$dens(xs)
    dados <- data.frame(x=xs,y=ys)
    dados <- dados[is.finite(dados$y) & dados$y>=0, ,drop=FALSE]
    validate(need(nrow(dados)>10,"A densidade não pôde ser desenhada."))
    selecionado <- switch(e$tipo,
      inferior = dados$x <= e$a,
      superior = dados$x >= e$a,
      intervalo = dados$x >= e$a & dados$x <= e$b,
      externa = dados$x <= e$a | dados$x >= e$b,
      quantil = dados$x <= e$a)
    dados$area <- ifelse(selecionado,dados$y,0)
    # Geom_ribbon aproxima a área visual. O valor impresso é EXATO via CDF.
    g <- ggplot(dados,aes(x=x)) +
      geom_ribbon(aes(ymin=0,ymax=area),fill=realce,alpha=.78) +
      geom_line(aes(y=y),linewidth=.85,color=azul) +
      geom_hline(yintercept=0,linewidth=.35,color="#718096") +
      labs(x="Valor da variável aleatória (x)",y="Densidade  f(x)") +
      theme_classic(base_size=13) +
      theme(plot.margin=margin(10,20,8,10),axis.title=element_text(color="#243750"),
            axis.text=element_text(color="#394B61")) +
      coord_cartesian(xlim=c(l,u),ylim=c(0,max(dados$y)*1.10),expand=FALSE)
    if (isTRUE(input$mostrar_media) && is.finite(m$media) && m$media >= l && m$media <= u)
      g <- g + geom_vline(xintercept=m$media,linetype="dotted",linewidth=.7,color="#CE7A30")
    if (isTRUE(input$mostrar_legenda)) {
      cortes <- if (e$tipo %in% c("intervalo","externa")) c(e$a,e$b) else e$a
      g <- g + geom_vline(xintercept=cortes,linetype="dashed",linewidth=.45,color="#237FA7")
    }
    if (input$tipo == "uniforme") {
      altura <- 1 / diff(m$suporte)
      g <- g + geom_segment(aes(x=m$suporte[1],xend=m$suporte[1],y=0,yend=altura),
                            inherit.aes=FALSE,color=azul,linewidth=.6) +
        geom_segment(aes(x=m$suporte[2],xend=m$suporte[2],y=0,yend=altura),
                     inherit.aes=FALSE,color=azul,linewidth=.6)
    }
    g
  })
  output$nota_grafico <- renderText({
    m <- modelo(); e <- evento()
    paste("A região azul representa uma aproximação gráfica da área selecionada;",
          "a probabilidade numérica é calculada pela função de distribuição acumulada.",
          "Nas caudas infinitas, a curva é exibida somente em uma janela que contém aproximadamente 99,9% da probabilidade.",
          m$observacao)
  })
  output$explicacao <- renderUI({
    m <- modelo()
    withMathJax(div(class="model-text",
      tags$p(tags$strong("Notação dos slides: "),math(m$notacao)),
      tags$p(tags$strong("Função densidade de probabilidade:")),
      tags$div(style="font-size:19px;overflow-x:auto;padding:8px 0;",math(m$formula)),
      tags$p(tags$strong("Domínio dos parâmetros e suporte: "),math(m$dominio)),
      tags$p(tags$strong("Significado dos parâmetros: "),m$parametros),
      tags$p(tags$strong("Esperança e variância: "),
             paste0("E(X) = ",fmt(m$media),"; Var(X) = ",fmt(m$variancia),".")),
      tags$p(tags$strong("Aplicações e interpretação: "),m$aplicacao),
      tags$p(class="note",m$observacao),
      tags$hr(),
      tags$p(class="note", "Para toda VAC, P(X = x₀) = 0 e P(a ≤ X ≤ b) = F(b) − F(a). A densidade f(x) pode exceder 1; sua integral total, e não sua altura máxima, deve ser igual a 1.")
    ))
  })
}

shinyApp(ui, server)
