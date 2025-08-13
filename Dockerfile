# Use the official Python 3.11 image as the base image
FROM python:3.11

# Set environment variables to non-interactive
ENV DEBIAN_FRONTEND=noninteractive

# Install Clang, sudo, and other necessary tools
RUN apt-get update && apt-get install -y \
    clang \
    sudo \
    wget \
    curl \
    git \
    python3-pip \
    zsh \
    autojump \
    && rm -rf /var/lib/apt/lists/*

# Create a new user 'user' with sudo access and zsh as the default shell
RUN useradd -m user -s $(which zsh) && echo "user:user" | chpasswd && adduser user sudo

# Switch to user
USER user
WORKDIR /home/user

# Copy the latest UV installer
COPY --from=ghcr.io/astral-sh/uv:0.8.9 /uv /uvx /bin/

# Create a Python virtual environment in /home/user/venv
RUN sh -c "wget -qO- https://astral.sh/uv/install.sh | sh"
RUN uv venv /home/user/venv -p python3.11
ENV PATH="/home/user/venv/bin:$PATH"

# Switch to root user
USER root

# Switch back to non-root user
USER user

# Set the working directory
WORKDIR /home/user

# Install Oh My Zsh
RUN sh -c "$(wget https://raw.github.com/ohmyzsh/ohmyzsh/master/tools/install.sh -O -)"

# Install zsh-autosuggestions plugin
RUN git clone https://github.com/zsh-users/zsh-autosuggestions ${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-autosuggestions

# Install zsh-syntax-highlighting plugin
RUN git clone https://github.com/zsh-users/zsh-syntax-highlighting.git ${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-syntax-highlighting

# Create a default .zshrc file with some suggested configurations
RUN echo 'export ZSH="$HOME/.oh-my-zsh"' >> $HOME/.zshrc && \
    echo 'ZSH_THEME="robbyrussell"' >> $HOME/.zshrc && \
    echo 'plugins=(git zsh-autosuggestions zsh-syntax-highlighting autojump)' >> $HOME/.zshrc && \
    echo 'source $ZSH/oh-my-zsh.sh' >> $HOME/.zshrc && \
    echo 'export HISTFILE=~/.zsh_history' >> $HOME/.zshrc && \
    echo 'export HISTSIZE=10000' >> $HOME/.zshrc && \
    echo 'export SAVEHIST=10000' >> $HOME/.zshrc && \
    echo 'setopt appendhistory' >> $HOME/.zshrc && \
    echo 'setopt histignorespace' >> $HOME/.zshrc && \
    echo 'setopt histignorealldups' >> $HOME/.zshrc && \
    echo 'source "$HOME/venv/bin/activate"' >> $HOME/.zshrc

# Set the default shell to zsh
SHELL ["/usr/bin/zsh", "-c"]

# Switch to root user
USER root

# Install Eigen3
RUN apt-get update && apt-get install -y libeigen3-dev

# Switch back to non-root user
USER user

# Install mdeq, accupy, cogent3 and jupyter / plotly stuff
RUN uv pip install "cogent3[extra]==2025.7.10a5"
RUN uv pip install "mdeq==2025.6.30"
RUN uv pip install accupy
# execute mdeq once to trigger cogent3 byte-compiling its things
RUN mdeq --version


# Setup editable install of the analysis code
WORKDIR /home/user/MutationDiseqAnalysis
ADD pyproject.toml README.md LICENSE ./
ADD src ./src
RUN uv pip install -e .

WORKDIR /home/user

# Start a terminal session using the zsh shell
CMD ["/usr/bin/zsh"]