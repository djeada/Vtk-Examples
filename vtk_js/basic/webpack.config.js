const path = require('path');

module.exports = {
  mode: 'development',
  entry: './src/main.js', // make sure this points to your main JavaScript file
  output: {
    path: path.resolve(__dirname, 'dist'),
    filename: 'bundle.js',
    libraryTarget: 'umd', // This exposes your bundle on different environments
    globalObject: 'this' // This ensures compatibility with both browser and Node environments
  },
  module: {
    rules: [
      {
        test: /\.js$/,
        exclude: /node_modules/,
        use: {
          loader: 'babel-loader',
          options: {
            presets: ['@babel/preset-env']
          }
        },
      },
    ],
  },
  resolve: {
    extensions: ['.js'], // Add '.ts' if you are using TypeScript
    fallback: {
      "stream": require.resolve("stream-browserify"),
      "buffer": require.resolve("buffer")
    }
  },
};

