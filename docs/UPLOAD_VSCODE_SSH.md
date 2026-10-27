# 用VS Code远程SSH终端上传已有GitHub仓库

适用于已经创建并有初始README的 `l1zhx/glide-fragment-screen`。本方法直接使用服务器上的Git，不依赖ChatGPT/Codex的GitHub应用授权。

## 1. 把发布包传到服务器

将 `GitHub_Glide_Workflow_20261027_v2.zip` 传到服务器 `~/github_publish/`。例如在VS Code连接服务器后，将本地ZIP拖入远程资源管理器中对应目录。后续命令全部在**服务器的远程终端**执行。

包中已经包含完整分析脚本和文档。不要把整个docking工作目录、配体库或轨迹复制进公开仓库。

## 2. 测试服务器到GitHub的认证

VS Code用SSH登录服务器，与服务器用SSH访问GitHub，是两条不同连接。先执行：

```bash
ssh -T git@github.com
```

成功消息应包含你的GitHub用户名l1zhx和successfully authenticated。GitHub不提供shell，成功测试也可能返回退出码1，不能仅按退出码判断失败。首次连接时核对显示的主机指纹。参见 [GitHub官方测试说明](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/testing-your-ssh-connection)。

若出现 `Permission denied (publickey)`，先按GitHub说明配置服务器的GitHub SSH密钥，完成后再继续。如果服务器已有可用的HTTPS凭据，也可将下面clone地址换成 `https://github.com/l1zhx/glide-fragment-screen.git`，继续使用原有认证方式。无需在脚本或仓库URL里写入token。

## 3. 克隆仓库，再复制发布文件

使用尚不存在的解压目录和克隆目录。仓库已经有初始提交，因此从clone开始，保留它的历史：

```bash
mkdir -p ~/github_publish
cd ~/github_publish

python -m zipfile -e GitHub_Glide_Workflow_20261027_v2.zip release_20261027
git clone git@github.com:l1zhx/glide-fragment-screen.git

cd glide-fragment-screen
git branch --show-current
git remote -v

cp -a ../release_20261027/glide-fragment-screen/. .
sha256sum -c MANIFEST.sha256
git status --short
```

分支应为main，origin应指向l1zhx/glide-fragment-screen，校验结果应全部OK。`cp`源目录末尾的 `/.` 会同时复制 `.github`、`.gitignore` 等隐藏文件；发布包不含 `.git`，所以不会覆盖克隆仓库的Git历史。它会用完整项目说明替换新仓库最初的简短README。

如果已经克隆过此仓库，直接进入那个干净的checkout并运行 `git pull --ff-only origin main` 后再复制；若有尚未提交的个人修改，先保留并处理这些修改，不要重复覆盖。

## 4. 检查后提交和推送

在克隆目录内：

```bash
git add README.md README.zh-CN.md LICENSE CHANGELOG.md requirements.txt VERSION MANIFEST.sha256 .gitignore .gitattributes .github scripts tests slurm config docs
git diff --cached --stat
```

检查暂存内容只包含项目代码、文档和合成示例图，再执行：

```bash
git commit -m "Add reproducible beta-catenin Glide screening workflow"
git push -u origin main
```

若提示缺少author identity，只在此仓库配置自己的姓名与GitHub邮箱后重试commit：

```bash
git config user.name "YOUR_GIT_NAME"
git config user.email "YOUR_GITHUB_EMAIL"
```

两个值需替换为自己的信息，可使用GitHub提供的noreply邮箱。如果push被拒绝为non-fast-forward，先检查远端新增提交并整合；不要使用force push来覆盖历史。这里不重新git init，也不使用unrelated-histories合并。

## 5. 查看发布结果

打开 [项目仓库](https://github.com/l1zhx/glide-fragment-screen)，确认根目录有README、scripts、docs和tests，并能显示示例图。进入Actions检查两种Python版本的公开测试任务；服务器docking不需要为了上传重新运行，公开CI也不会运行商业软件。

## 6. 将项目显示在个人主页

GitHub个人主页README使用账号同名的公开仓库 `l1zhx/l1zhx`；它与Glide源码仓库分开。可以先在个人主页Pinned中固定Glide仓库。若要使用本包的完整主页，创建同名仓库并勾选添加README后：

```bash
cd ~/github_publish
git clone git@github.com:l1zhx/l1zhx.git profile_checkout
cd profile_checkout
cat README.md
```

如果只有初始化标题，复制交付的完整主页后提交：

```bash
cp ../release_20261027/l1zhx/README.md README.md
git diff -- README.md
git add README.md
git commit -m "Introduce molecular screening projects"
git push -u origin main
```

如果主页已经有自己的介绍，用VS Code将 [profile_snippet.md](profile_snippet.md) 的项目段落合并到现有README，再提交；不要直接覆盖。发布后的主页会同时链接AF3与Glide两个项目。
