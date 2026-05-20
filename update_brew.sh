#!/bin/bash

set -e

if [ -z "$1" ]; then
  echo "Usage: ./update_brew.sh <version>"
  echo "Example: ./update_brew.sh 1.0.1"
  exit 1
fi

VERSION="$1"
TAG="v$VERSION"

PROJECT_DIR="$HOME/taykit"
TAP_DIR="$HOME/homebrew-taykit"

TARBALL_NAME="taykit-macos-arm64.tar.gz"
RELEASE_URL="https://github.com/richardmarktaylor/taykit/releases/download/$TAG/$TARBALL_NAME"

cd "$PROJECT_DIR"

echo "Building taykit..."
./build.sh

echo "Creating tarball..."
cd dist
rm -f "$TARBALL_NAME"
tar -czf "$TARBALL_NAME" taykit

SHA256=$(shasum -a 256 "$TARBALL_NAME" | awk '{print $1}')

echo "SHA256: $SHA256"

cd "$PROJECT_DIR"

echo "Committing project changes..."
git add .
git commit -m "Release $TAG" || echo "No project changes to commit"
git push

echo "Creating GitHub release..."
gh release create "$TAG" \
  "dist/$TARBALL_NAME" \
  --title "taykit $TAG" \
  --notes "Release $TAG" || echo "Release may already exist"

echo "Updating Homebrew formula..."

cd "$TAP_DIR"

cat > Formula/taykit.rb <<EOF
class Taykit < Formula
  desc "Taylor bioinformatics command-line toolkit"
  homepage "https://github.com/richardmarktaylor/taykit"
  url "$RELEASE_URL"
  sha256 "$SHA256"
  version "$VERSION"

  def install
    bin.install "taykit"
  end

  test do
    system "#{bin}/taykit", "--help"
  end
end
EOF

git add Formula/taykit.rb
git commit -m "Update taykit to $VERSION"
git push

echo
echo "Done."
echo
echo "Users can now run:"
echo "  brew update"
echo "  brew upgrade taykit"
