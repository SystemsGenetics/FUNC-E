### Image for FUNC-E perl
FROM rocker/rstudio
LABEL maintainer="Reed Bender <mrbende@clemson.edu>"

# Download dependencies for FUNC_E
RUN apt-get update && \
cpan App::cpanminus && \
cpan inc::latest && \
cpan IPC::Run \
cpan Getopt::Long && \
cpan Text::NSP::Measures::2D::Fisher::right && \
cpan List::Util && \
cpan Math::Complex && \
cpan Math::BigFloat && \
cpan Statistics::R 
